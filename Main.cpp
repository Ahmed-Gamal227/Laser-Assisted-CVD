#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <random>
#include <limits>
#include <algorithm>

// Constants
const double kB = 8.617333262145e-5; // Boltzmann constant in eV/K
const double NA = 6.02214076e23;     // Avogadro's number

// Simulation parameters
const int GRID_SIZE = 100;           // Size of the simulation grid
const double CELL_SIZE = 1e-6;       // 1 micron cell size in meters
const double SIM_TIME = 1.0;         // Total simulation time in seconds
const double DT = 0.01;              // Time step in seconds
const double OUTPUT_INTERVAL = 0.1;  // Time between outputs

// Material properties
const double SUBSTRATE_THERMAL_CONDUCTIVITY = 150.0; // W/mK (e.g., silicon)
const double SPECIFIC_HEAT_CAPACITY = 700.0;         // J/kgK
const double MATERIAL_DENSITY = 2330.0;              // kg/m^3 (silicon)
const double INITIAL_TEMP = 300.0;                   // Initial temperature (K)

// LCVD parameters
const double LASER_POWER = 1.0;      // Laser power in Watts
const double LASER_SPOT_RADIUS = 5e-6; // Laser spot radius in meters (5 microns)
const double ACTIVATION_ENERGY = 1.2; // Activation energy in eV
const double PRE_EXPONENTIAL_FACTOR = 1e12; // Pre-exponential factor in 1/s
const double GAS_PRESSURE = 100.0;   // Gas pressure in Pa
const double STICKING_COEFFICIENT = 0.1; // Probability of gas molecule sticking
const double MIN_DEPOSITION_RATE = 1e-20; // Minimum deposition rate

class LCVD_Simulation {
private:
    std::vector<std::vector<double>> temperature;
    std::vector<std::vector<double>> thickness;
    std::vector<std::vector<double>> gas_concentration;

    double time;
    std::mt19937 gen;
    std::uniform_real_distribution<> dis;

public:
    LCVD_Simulation() : time(0.0), gen(std::random_device{}()), dis(0.0, 1.0) {
        // Initialize grids with proper dimensions
        temperature.resize(GRID_SIZE, std::vector<double>(GRID_SIZE, INITIAL_TEMP));
        thickness.resize(GRID_SIZE, std::vector<double>(GRID_SIZE, 0.0));
        gas_concentration.resize(GRID_SIZE, std::vector<double>(GRID_SIZE, 1.0));
    }

    // Gaussian laser beam profile with bounds checking
    double laser_heat_source(int x, int y) {
        if (x < 0 || x >= GRID_SIZE || y < 0 || y >= GRID_SIZE) return 0.0;

        double x_pos = (x - GRID_SIZE / 2.0) * CELL_SIZE;
        double y_pos = (y - GRID_SIZE / 2.0) * CELL_SIZE;
        double r2 = x_pos * x_pos + y_pos * y_pos;

        double intensity = (LASER_POWER / (3.14 * LASER_SPOT_RADIUS * LASER_SPOT_RADIUS)) *
            exp(-r2 / (LASER_SPOT_RADIUS * LASER_SPOT_RADIUS));

        return std::isfinite(intensity) ? intensity : 0.0;
    }

    // Safe Arrhenius reaction rate calculation
    double deposition_rate(double T) {
        if (T <= 300.0) return 0.0;

        try {
            double rate = PRE_EXPONENTIAL_FACTOR * exp(-ACTIVATION_ENERGY / (kB * T));
            return std::isfinite(rate) ? std::max(rate, MIN_DEPOSITION_RATE) : 0.0;
        }
        catch (...) {
            return 0.0;
        }
    }

    // Update temperature field with stability checks
    void update_temperature() {
        std::vector<std::vector<double>> new_temp = temperature;

        double alpha = SUBSTRATE_THERMAL_CONDUCTIVITY / (SPECIFIC_HEAT_CAPACITY * MATERIAL_DENSITY);
        double dt_over_dx2 = alpha * DT / (CELL_SIZE * CELL_SIZE);

        // Ensure stability condition (dt <= dx^2/(4*alpha))
        if (dt_over_dx2 > 0.25) {
            std::cerr << "Warning: Unstable temperature update, reducing time step\n";
            dt_over_dx2 = 0.25;
        }

        for (int i = 1; i < GRID_SIZE - 1; ++i) {
            for (int j = 1; j < GRID_SIZE - 1; ++j) {
                double heat_source = laser_heat_source(i, j);

                double temp_update = temperature[i][j] + dt_over_dx2 * (
                    temperature[i + 1][j] + temperature[i - 1][j] +
                    temperature[i][j + 1] + temperature[i][j - 1] -
                    4.0 * temperature[i][j]) +
                    (heat_source * DT) / (SPECIFIC_HEAT_CAPACITY * MATERIAL_DENSITY * CELL_SIZE);

                // Validate temperature
                if (!std::isfinite(temp_update) || temp_update < 0) {
                    new_temp[i][j] = INITIAL_TEMP;
                }
                else {
                    new_temp[i][j] = temp_update;
                }
            }
        }

        // Boundary conditions
        for (int i = 0; i < GRID_SIZE; ++i) {
            new_temp[i][0] = new_temp[i][GRID_SIZE - 1] = INITIAL_TEMP;
            new_temp[0][i] = new_temp[GRID_SIZE - 1][i] = INITIAL_TEMP;
        }

        temperature = new_temp;
    }

    // Update deposition thickness with safety checks
    void update_deposition() {
        for (int i = 0; i < GRID_SIZE; ++i) {
            for (int j = 0; j < GRID_SIZE; ++j) {
                double rate = deposition_rate(temperature[i][j]);
                double prob = STICKING_COEFFICIENT * gas_concentration[i][j];

                if (dis(gen) < prob * rate * DT) {
                    double deposit = 1e-10; // 0.1 nm per deposition event
                    if (std::isfinite(deposit) && deposit > 0) {
                        thickness[i][j] += deposit;
                    }
                }

                // Ensure thickness doesn't become negative
                thickness[i][j] = std::max(0.0, thickness[i][j]);
            }
        }
    }

    // Update gas concentration with diffusion stability
    void update_gas_concentration() {
        std::vector<std::vector<double>> new_conc = gas_concentration;
        double D = 1e-5; // Diffusion coefficient (m^2/s)
        double dt_over_dx2 = D * DT / (CELL_SIZE * CELL_SIZE);

        // Stability check
        if (dt_over_dx2 > 0.25) {
            std::cerr << "Warning: Unstable gas diffusion, reducing time step\n";
            dt_over_dx2 = 0.25;
        }

        for (int i = 1; i < GRID_SIZE - 1; ++i) {
            for (int j = 1; j < GRID_SIZE - 1; ++j) {
                double conc_update = gas_concentration[i][j] + dt_over_dx2 * (
                    gas_concentration[i + 1][j] + gas_concentration[i - 1][j] +
                    gas_concentration[i][j + 1] + gas_concentration[i][j - 1] -
                    4.0 * gas_concentration[i][j]);

                // Depletion due to deposition
                double dep_rate = deposition_rate(temperature[i][j]);
                conc_update -= dep_rate * DT * STICKING_COEFFICIENT;

                // Validate concentration
                if (!std::isfinite(conc_update)) {
                    new_conc[i][j] = 1.0;
                }
                else {
                    new_conc[i][j] = std::max(0.0, std::min(1.0, conc_update));
                }
            }
        }

        // Boundary conditions
        for (int i = 0; i < GRID_SIZE; ++i) {
            new_conc[i][0] = new_conc[i][GRID_SIZE - 1] = 1.0;
            new_conc[0][i] = new_conc[GRID_SIZE - 1][i] = 1.0;
        }

        gas_concentration = new_conc;
    }

    void run_simulation() {
        std::ofstream temp_file("temperature.csv");
        std::ofstream thickness_file("thickness.csv");

        if (!temp_file.is_open() || !thickness_file.is_open()) {
            std::cerr << "Error opening output files!\n";
            return;
        }

        int output_count = 0;

        for (time = 0; time <= SIM_TIME; time += DT) {
            update_temperature();
            update_gas_concentration();
            update_deposition();

            // Output data periodically
            if (fmod(time, OUTPUT_INTERVAL) < DT / 2.0) {
                save_slice(temp_file, temperature);
                save_slice(thickness_file, thickness);
                output_count++;

                std::cout << "Progress: " << (time / SIM_TIME) * 100 << "% \r";
                std::cout.flush();
            }
        }

        temp_file.close();
        thickness_file.close();
        std::cout << "\nSimulation completed. " << output_count << " snapshots saved.\n";
    }

    void save_slice(std::ofstream& file, const std::vector<std::vector<double>>& data) {
        for (const auto& row : data) {
            for (size_t j = 0; j < row.size(); ++j) {
                file << row[j];
                if (j < row.size() - 1) file << ",";
            }
            file << "\n";
        }
        file << "\n";
    }
};

int main() {
    try {
        LCVD_Simulation sim;
        sim.run_simulation();
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
} 