#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

using namespace std;

int main() {
  // --- PARAMETERS ---
  const int nx = 256, ny = 256;
  const int pad_x = nx + 2,
            pad_y = ny + 2; // Including buffer layers (x1-1 to x2+1)

  const int np_initial = 14;
  const int ssteps = 10;

  const double dh = 1.0;
  const double gamma = 1.5;
  const double W = 1.0;
  const double eps_sq = 1.0;
  const double sLdt = 0.1;

  // Boundary Condition Parameters
  const double cx = nx / 2.0;
  const double cy = ny / 2.0;
  const double R = (nx / 2.0) * 0.85;
  const double R_sq = R * R;

  auto idx2 = [pad_x](int x, int y) { return (y * pad_x) + x; };

  auto idx3 = [pad_x](int x, int y, int p, int num_p) {
    return ((y * pad_x) + x) * num_p + p;
  };

  // --- INITIAL VORONOI TESSELLATION ---
  cout << "Initializing Voronoi tessellation..." << endl;
  vector<int> featureID(pad_x * pad_y, 0);
  vector<vector<double>> centers(np_initial, vector<double>(2));
  /*
  mt19937 gen(44392); // Replicating the Fortran random seed
  uniform_real_distribution<double> dist(0.0, 1.0);

  double length[2] = {nx * dh, ny * dh};
  for (int ip = 0; ip < np_initial; ++ip) {
    centers[ip][0] = dist(gen) * length[0];
    centers[ip][1] = dist(gen) * length[1];
  }
  */

  // Fixed 4 grain case, ensure np = 8
  // Sort this out later, currently a little odd
  /*
  centers[0][0] = 0.4 * nx; // 1x
  centers[0][1] = 0.4 * ny; // 1y
  centers[1][0] = 0.6 * nx; // 2x
  centers[1][1] = 0.6 * ny; // 2y
  centers[2][0] = 0.4 * nx; // 3x
  centers[2][1] = 0.6 * ny; // 3y
  centers[3][0] = 0.6 * nx; // 4x
  centers[3][1] = 0.4 * ny; // 4y
  centers[4][0] = 0.1 * nx; // 1x
  centers[4][1] = 0.1 * ny; // 1y
  centers[5][0] = 0.9 * nx; // 2x
  centers[5][1] = 0.9 * ny; // 2y
  centers[6][0] = 0.1 * nx; // 3x
  centers[6][1] = 0.9 * ny; // 3y
  centers[7][0] = 0.9 * nx; // 4x
  centers[7][1] = 0.1 * ny; // 4y

  // Bicrystal case, np = 14

  centers[0][0] = 0.5 * nx;   // 1x
  centers[0][1] = 0.475 * ny; // 1y
  centers[1][0] = 0.5 * nx;
  centers[1][1] = 0.525 * ny;

  centers[2][0] = 0.5 * nx;
  centers[2][1] = 0.05 * ny;
  centers[3][0] = 0.05 * nx;
  centers[3][1] = 0.5 * ny;
  centers[4][0] = 0.5 * nx;
  centers[4][1] = 0.95 * ny;
  centers[5][0] = 0.95 * nx;
  centers[5][1] = 0.5 * ny;

  centers[6][0] = 0.818198 * nx;
  centers[6][1] = 0.818198 * ny;
  centers[7][0] = 0.818198 * nx;
  centers[7][1] = 0.181802 * ny;
  centers[8][0] = 0.181802 * nx;
  centers[8][1] = 0.181802 * ny;
  centers[9][0] = 0.181802 * nx;
  centers[9][1] = 0.818198 * ny;

  centers[10][0] = 0.9 * nx;
  centers[10][1] = 0.9 * ny;
  centers[11][0] = 0.9 * nx;
  centers[11][1] = 0.1 * ny;
  centers[12][0] = 0.1 * nx;
  centers[12][1] = 0.1 * ny;
  centers[13][0] = 0.1 * nx;
  centers[13][1] = 0.9 * ny;
  */

  for (int iy = 0; iy < pad_y; ++iy) {
    for (int ix = 0; ix < pad_x; ++ix) {
      double dist2min = 1e30; // huge value
      double coord[2] = {ix * dh, iy * dh};

      for (int ip = 0; ip < np_initial; ++ip) {
        double dist2 = 0.0;
        for (int di = 0; di < 2; ++di) {
          double d = coord[di] - centers[ip][di];
          // Periodic distance (commented out in Fortran source, leaving
          // simple distance here)
          dist2 += d * d;
        }
        if (dist2 < dist2min) {
          dist2min = dist2;
          featureID[idx2(ix, iy)] = ip; // 0-indexed ip
        }
      }
    }
  }

  // --- KNOCK OUT GRAINS ---
  cout << "Knocking out invalid grains..." << endl;
  vector<int> index_to_featureID;
  int npeff = 0;

  for (int ip = 0; ip < np_initial; ++ip) {
    int count_mask = 0;
    for (int v : featureID) {
      if (v == ip)
        count_mask++;
    }

    double dist2_center = pow(centers[ip][0] - (nx / 2.0), 2) +
                          pow(centers[ip][1] - (nx / 2.0), 2);
    double radius_limit = pow((nx / 2.0) * 0.85, 2);

    if (count_mask < 300 || dist2_center > radius_limit) {
      for (int &v : featureID) {
        if (v == ip)
          v = np_initial; // mapped to invalid grain
      }
    } else {
      index_to_featureID.push_back(ip);
      npeff++;
    }
  }

  index_to_featureID.push_back(
      np_initial); // makes it a meatball, remove for full structure
  npeff++;

  cout << "Grains reduced from " << np_initial << " to " << npeff << endl;

  // Fixed case where there are two grains inside a circle
  centers[0][0] = 0.4 * nx;
  centers[0][1] = 0.5 * ny;
  centers[1][0] = 0.6 * nx;
  centers[1][1] = 0.5 * ny;

  // --- INITIAL VORONOI TESSELATION --- //
  for (int iy = 0; iy < pad_y; ++iy) {
    for (int ix = 0; ix < pad_x; ++ix) {
      double dist2_center =
          (ix * dh - cx) * (ix * dh - cx) + (iy * dh - cy) * (iy * dh - cy);

      // If outside radial limit, assign to boundary mask phase (essentially the
      // shadow realm)
      if (dist2_center > R_sq) {
        featureID[idx2(ix, iy)] = np_initial; // 0:np_initial-1 are grains,
                                              // np_initial is the shadow realm
      } else {
        // Standard Voronoi Algorithm
        double dist2min = 1e30;
        for (int ip = 0; ip < np_initial; ++ip) {
          double dist2 = 0.0;
          double dx = (ix * dh) - centers[ip][0];
          double dy = (iy * dh) - centers[ip][1];
          dist2 = dx * dx + dy * dy;

          if (dist2 < dist2min) {
            dist2min = dist2;
            featureID[idx2(ix, iy)] = ip;
          }
        }
      }
    }
  }

  // --- INITIALIZE ORDER PARAMETERS (PHI) --- //
  /*
  vector<double> Phi(pad_x * pad_y * npeff, 0.0);
  for (int ip = 0; ip < npeff; ++ip) {
    for (int iy = 0; iy < pad_y; ++iy) {
      for (int ix = 0; ix < pad_x; ++ix) {
        if (featureID[idx2(ix, iy)] == index_to_featureID[ip]) {
          Phi[idx3(ix, iy, ip, npeff)] = 1.0;
        }
      }
    }
  }

  featureID.clear(); // Free memory
  index_to_featureID.clear();

  vector<double> Phi_new = Phi;
  vector<double> smsq(pad_x * pad_y, 0.0);
  */

  // --- INITIALIZE ORDER PARAMETERS (PHI) --- //
  int npeff = np_initial + 1;
  cout << "Total effective phases (grains + boundary): " << npeff << endl;

  vector<int> index_to_featureID(npeff);
  iota(index_to_featureID.begin(), index_to_featureID.end(),
       0); // Maps 0, 1, 2 smoothly

  vector<double> Phi(pad_x *)

      // --- TIME EVOLUTION ---
      cout
      << "Starting time integration (" << ssteps << " steps)..." << endl;
  for (int it = 0; it <= ssteps; ++it) {

    // Compute smsq = sum(Phi^2)
    fill(smsq.begin(), smsq.end(), 0.0);
    for (int iy = 1; iy <= ny; ++iy) {
      for (int ix = 1; ix <= nx; ++ix) {
        double sum_sq = 0.0;
        for (int ip = 0; ip < npeff; ++ip) {
          double p = Phi[idx3(ix, iy, ip, npeff)];
          sum_sq += p * p;
        }
        smsq[idx2(ix, iy)] = sum_sq;
      }
    }

    // Evolve Phi
    for (int iy = 1; iy <= ny; ++iy) {
      for (int ix = 1; ix <= nx; ++ix) {
        for (int ip = 0; ip < npeff; ++ip) {
          double p = Phi[idx3(ix, iy, ip, npeff)];
          double dfdphi =
              p * p * p - p + 2.0 * gamma * p * (smsq[idx2(ix, iy)] - p * p);

          double lapl = (Phi[idx3(ix + 1, iy, ip, npeff)] +
                         Phi[idx3(ix - 1, iy, ip, npeff)] +
                         Phi[idx3(ix, iy + 1, ip, npeff)] +
                         Phi[idx3(ix, iy - 1, ip, npeff)] - 4.0 * p) /
                        (dh * dh);

          Phi_new[idx3(ix, iy, ip, npeff)] =
              p - sLdt * (W * dfdphi - eps_sq * lapl);
        }
      }
    }

    // Apply Periodic BCs in buffer layers
    for (int ip = 0; ip < npeff; ++ip) {
      for (int iy = 0; iy < pad_y; ++iy) {
        Phi_new[idx3(0, iy, ip, npeff)] = Phi_new[idx3(nx, iy, ip, npeff)];
        Phi_new[idx3(nx + 1, iy, ip, npeff)] = Phi_new[idx3(1, iy, ip, npeff)];
      }
      for (int ix = 0; ix < pad_x; ++ix) {
        Phi_new[idx3(ix, 0, ip, npeff)] = Phi_new[idx3(ix, ny, ip, npeff)];
        Phi_new[idx3(ix, ny + 1, ip, npeff)] = Phi_new[idx3(ix, 1, ip, npeff)];
      }
    }
    Phi = Phi_new;

    // --- DYNAMIC GRAIN PRUNING ---
    if (it % 100 == 0 || it == ssteps) {
      vector<int> index_to_index;
      for (int ip = 0; ip < npeff; ++ip) {
        double vol_sum = 0.0;
        for (int iy = 0; iy < pad_y; ++iy) {
          for (int ix = 0; ix < pad_x; ++ix) {
            vol_sum += Phi[idx3(ix, iy, ip, npeff)];
          }
        }
        if (vol_sum > 10.0)
          index_to_index.push_back(ip);
      }

      int npeff_fin = index_to_index.size();
      if (npeff_fin < npeff) {
        cout << "Step " << it << ": Reducing grains from " << npeff << " to "
             << npeff_fin << endl;
        vector<double> Phi_shrink(pad_x * pad_y * npeff_fin, 0.0);

        for (int ip = 0; ip < npeff_fin; ++ip) {
          int old_ip = index_to_index[ip];
          for (int iy = 0; iy < pad_y; ++iy) {
            for (int ix = 0; ix < pad_x; ++ix) {
              Phi_shrink[idx3(ix, iy, ip, npeff_fin)] =
                  Phi[idx3(ix, iy, old_ip, npeff)];
            }
          }
        }
        Phi = Phi_shrink;
        Phi_new = Phi;
        npeff = npeff_fin;
      } else if (it % 500 == 0) {
        cout << "Step " << it << " completed." << endl;
      }
    }
  }

  // --- POST-PROCESSING ---
  cout << "Normalizing and writing outputs..." << endl;
  for (double &val : Phi) {
    if (val < 0.0)
      val = 0.0;
    if (val > 1.0)
      val = 1.0;
  }

  for (int iy = 0; iy < pad_y; ++iy) {
    for (int ix = 0; ix < pad_x; ++ix) {
      double local_sum = 0.0;
      for (int ip = 0; ip < npeff; ++ip)
        local_sum += Phi[idx3(ix, iy, ip, npeff)];
      if (local_sum > 0) {
        for (int ip = 0; ip < npeff; ++ip)
          Phi[idx3(ix, iy, ip, npeff)] /= local_sum;
      }
    }
  }

  // Recalculate smsq (Cross terms only for grain boundaries)
  fill(smsq.begin(), smsq.end(), 0.0);
  for (int iy = 1; iy <= ny; ++iy) {
    for (int ix = 1; ix <= nx; ++ix) {
      double sum_cross = 0.0;
      for (int ip = 0; ip < npeff; ++ip) {
        for (int jp = 0; jp < npeff; ++jp) {
          if (ip != jp) {
            sum_cross +=
                Phi[idx3(ix, iy, ip, npeff)] * Phi[idx3(ix, iy, jp, npeff)];
          }
        }
      }
      smsq[idx2(ix, iy)] = sum_cross;
    }
  }

  // --- FILE EXPORT ---
  // Extract interior grid for outputs (ignoring pad layers, mapping to exact
  // dimensions make_props expects)

  // 1. grain_vis.dat
  vector<double> smsq_out(nx * ny);
  for (int iy = 1; iy <= ny; ++iy)
    for (int ix = 1; ix <= nx; ++ix)
      smsq_out[(iy - 1) * nx + (ix - 1)] = smsq[idx2(ix, iy)];

  ofstream outVis("grain_vis.dat", ios::binary);
  if (outVis.is_open()) {
    outVis.write(reinterpret_cast<char *>(smsq_out.data()),
                 smsq_out.size() * sizeof(double));
    outVis.close();
  }

  // 2. grains.dat (This must perfectly match the 4D interior space for your
  // other C++ script)
  vector<double> Phi_out(nx * ny * npeff);

  for (int iy = 1; iy <= ny; ++iy) {
    for (int ix = 1; ix <= nx; ++ix) {
      for (int ip = 0; ip < npeff; ++ip) {
        int out_idx = ((ix - 1) * ny + (iy - 1)) * npeff +
                      ip; // matches idx4D in make_props
        Phi_out[out_idx] = Phi[idx3(ix, iy, ip, npeff)];
      }
    }
  }

  ofstream outGrains("grains.dat", ios::binary);
  if (outGrains.is_open()) {
    outGrains.write(reinterpret_cast<char *>(Phi_out.data()),
                    Phi_out.size() * sizeof(double));
    outGrains.close();
  }

  cout << "Finished successfully!" << endl;
  return 0;
}