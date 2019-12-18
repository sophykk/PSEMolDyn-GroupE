
#include "FileReader.h"
#include "outputWriter/XYZWriter.h"
#include "utils/ArrayUtils.h"

#include <iostream>
#include <list>

using namespace std;

/**** forward declaration of the calculation functions ****/

/**
 * calculate the force for all particles
 */
void calculateF();

/**
 * calculate the position for all particles
 */
void calculateX();

/**
 * calculate the position for all particles
 */
void calculateV();

/**
 * plot the particles to a xyz-file
 */
void plotParticles(int iteration);

constexpr double start_time = 0;
constexpr double end_time = 1000;
constexpr double delta_t = 0.014;

std::list<Particle> particles;

int main(int argc, char *argsv[]) {

  cout << "Hello from MolSim for PSE!" << endl;
  if (argc != 2) {
    cout << "Errounous programme call! " << endl;
    cout << "./molsym filename" << endl;
  }

  FileReader fileReader;
  fileReader.readFile(particles, argsv[1]);

  double current_time = start_time;

  int iteration = 0;

  // for this loop, we assume: current x, current f and current v are known
  while (current_time < end_time) {
    // calculate new x
    calculateX();

    // calculate new f
    calculateF();
    // calculate new v
    calculateV();

    iteration++;
    if (iteration % 10 == 0) {
      plotParticles(iteration);
    }
    cout << "Iteration " << iteration << " finished." << endl;

    current_time += delta_t;
  }

  cout << "output written. Terminating..." << endl;
  return 0;
}

void calculateF() {
  list<Particle>::iterator iterator;
  iterator = particles.begin();

  for (auto &p1 : particles) {
    for (auto &p2 : particles) {
      // @TODO: insert calculation of force here!
    }
  }
}

void calculateX() {
  for (auto &p : particles) {
      // @TODO: insert calculation of force here!
  }
}

void calculateV() {
  for (auto &p : particles) {
    // @TODO: insert calculation of force here!
  }
}

void plotParticles(int iteration) {

  string out_name("MD_vtk");

  outputWriter::XYZWriter writer;
  writer.plotParticles(particles, out_name, iteration);
}
