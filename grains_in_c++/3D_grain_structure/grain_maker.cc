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
  const int nx = 97, ny = 97, nz = 97;
  const int pad_x = nx + 2, pad_y = ny + 2,
            pad_z = nz + 2; // Including buffer layers (x1-1 to x2+1)
  const int extrax = 16;
  const int innerx = nx - 2*extrax, innery = ny - 2*extrax, innerz = nz - 2*extrax; 
  const int np_initial = 200;
  const int ssteps = 700;

  const double dh = 1.0;
  const double gamma = 1.5;
  const double W = 1.0;
  const double eps_sq = 1.0;
  const double sLdt = 0.08;

  // Helper lambdas for multi-dimensional 1D indexing
  // Layout: [z][y][x][p] to optimize cache locality during the spatial stencils
  auto idx3 = [pad_x, pad_y](int x, int y, int z) {
    return (z * pad_y * pad_x) + (y * pad_x) + x;
  };

  auto idx4 = [pad_x, pad_y](int x, int y, int z, int p, int num_p) {
    return ((z * pad_y * pad_x) + (y * pad_x) + x) * num_p + p;
  };
  
  auto idx3_small = [innerx, innery](int x, int y, int z) {
    return (z * innery * innerx) + (y * innerx) + x;
  };
  
  auto idx4_small = [innerx, innery](int x, int y, int z, int p, int num_p) {
    return ((z * innery * innerx) + (y * innerx) + x) * num_p + p;
  };

  // --- INITIAL VORONOI TESSELLATION ---
  cout << "Initializing Voronoi tessellation..." << endl;
  vector<int> featureID(pad_x * pad_y * pad_z, 0);
  vector<vector<double>> centers(np_initial, vector<double>(3));

  mt19937 gen(44392); // Replicating the Fortran random seed
  uniform_real_distribution<double> dist(0.0, 1.0);

  double length[3] = {nx * dh, ny * dh, nz * dh};
  for (int ip = 0; ip < np_initial; ++ip) {
    centers[ip][0] = dist(gen) * length[0];
    centers[ip][1] = dist(gen) * length[1];
    centers[ip][2] = dist(gen) * length[2];
  }
  
  // ip = 0 is now the big grain
  centers[0][0] = 0.0;
  centers[0][1] = 0.0;
  centers[0][2] = 0.0;

  for (int iz = 0; iz < pad_z; ++iz) {
    for (int iy = 0; iy < pad_y; ++iy) {
      for (int ix = 0; ix < pad_x; ++ix) {
        double dist2min = 1e30; // huge value
        double coord[3] = {ix * dh, iy * dh, iz * dh};

        for (int ip = 0; ip < np_initial; ++ip) {
          double dist2 = 0.0;
          for (int di = 0; di < 3; ++di) {
            double d = coord[di] - centers[ip][di];
            // Periodic distance (commented out in Fortran source, leaving
            // simple distance here)
            dist2 += d * d;
          }
          if (dist2 < dist2min) {
            dist2min = dist2;
            featureID[idx3(ix, iy, iz)] = ip; // 0-indexed ip
          }
        }
      }
    }
  }
  
  // --- KNOCK OUT GRAINS ---
  cout << "Knocking out invalid grains..." << endl;
  vector<int> index_to_featureID;
  int npeff = 0;

  for (int ip = 1; ip < np_initial; ++ip) {
    int count_mask = 0;
    for (int v : featureID) {
      if (v == ip)
        count_mask++;
    }

    if (count_mask < 100) {
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

  // --- INITIALIZE ORDER PARAMETERS (PHI) ---
  vector<double> Phi(pad_x * pad_y * pad_z * npeff, 0.0);
  for (int iz = 0; iz < pad_z; ++iz) {
    for (int iy = 0; iy < pad_y; ++iy) {
      for (int ix = 0; ix < pad_x; ++ix) {
        double dist = sqrt(pow(ix - (nx / 2.0), 2) +
                                 pow(iy - (ny / 2.0), 2) +
                                 pow(iz - (nz / 2.0), 2));
        double distfn = 0.5*(1.0+tanh((innerx/2*1.2-dist)/2.0));
        Phi[idx4(ix, iy, iz, 0, npeff)] = 1.0-distfn;
        for (int ip = 0; ip < npeff; ++ip) {
          if (featureID[idx3(ix, iy, iz)] == ip) {
            Phi[idx4(ix, iy, iz, ip, npeff)] = distfn;
          }
        }
      }
    }
  }

  featureID.clear(); // Free memory
  index_to_featureID.clear();
  
  // Applying Dirichlet BCs
  for (int ip = 0; ip < npeff; ++ip) {
    for (int iz = 0; iz < pad_z; ++iz) {
      for (int iy = 0; iy < pad_y; ++iy) {
        Phi[idx4(0, iy, iz, ip, npeff)] = 0.0;
        Phi[idx4(nx + 1, iy, iz, ip, npeff)] = 0.0;
      }
      for (int ix = 0; ix < pad_x; ++ix) {
        Phi[idx4(ix, 0, iz, ip, npeff)] = 0.0;
        Phi[idx4(ix, ny + 1, iz, ip, npeff)] = 0.0;
      }
    }
    for (int iy = 0; iy < pad_y; ++iy) {
      for (int ix = 0; ix < pad_x; ++ix) {
        Phi[idx4(ix, iy, 0, ip, npeff)] = 0.0;
        Phi[idx4(ix, iy, nz + 1, ip, npeff)] = 0.0;
      }
    }
  }

  vector<double> Phi_new = Phi;
  vector<double> smsq(pad_x * pad_y * pad_z, 0.0);

  

  // --- TIME EVOLUTION ---
  cout << "Starting time integration (" << ssteps << " steps)..." << endl;
  for (int it = 0; it <= ssteps; ++it) {

    // Compute smsq = sum(Phi^2)
    fill(smsq.begin(), smsq.end(), 0.0);
    for (int iz = 1; iz <= nz; ++iz) {
      for (int iy = 1; iy <= ny; ++iy) {
        for (int ix = 1; ix <= nx; ++ix) {
          double sum_sq = 0.0;
          for (int ip = 0; ip < npeff; ++ip) {
            double p = Phi[idx4(ix, iy, iz, ip, npeff)];
            sum_sq += p * p;
          }
          smsq[idx3(ix, iy, iz)] = sum_sq;
        }
      }
    }

    // Evolve Phi
    for (int iz = 1; iz <= nz; ++iz) {
      for (int iy = 1; iy <= ny; ++iy) {
        for (int ix = 1; ix <= nx; ++ix) {
          for (int ip = 0; ip < npeff; ++ip) {
            double p = Phi[idx4(ix, iy, iz, ip, npeff)];
            double dfdphi = p * p * p - p +
                            2.0 * gamma * p * (smsq[idx3(ix, iy, iz)] - p * p);

            double lapl = (Phi[idx4(ix + 1, iy, iz, ip, npeff)] +
                           Phi[idx4(ix - 1, iy, iz, ip, npeff)] +
                           Phi[idx4(ix, iy + 1, iz, ip, npeff)] +
                           Phi[idx4(ix, iy - 1, iz, ip, npeff)] +
                           Phi[idx4(ix, iy, iz + 1, ip, npeff)] +
                           Phi[idx4(ix, iy, iz - 1, ip, npeff)] - 6.0 * p) /
                          (dh * dh);

            Phi_new[idx4(ix, iy, iz, ip, npeff)] =
                p - sLdt * (W * dfdphi - eps_sq * lapl);
          }
        }
      }
    }

    Phi = Phi_new;

    // --- DYNAMIC GRAIN PRUNING ---
    if (it % 100 == 0 || it == ssteps) {
      vector<int> index_to_index;
      for (int ip = 0; ip < npeff; ++ip) {
        double vol_sum = 0.0;
        for (int iz = 0; iz < pad_z; ++iz) {
          for (int iy = 0; iy < pad_y; ++iy) {
            for (int ix = 0; ix < pad_x; ++ix) {
              vol_sum += Phi[idx4(ix, iy, iz, ip, npeff)];
            }
          }
        }
        if (vol_sum > 20.0)
          index_to_index.push_back(ip);
      }

      int npeff_fin = index_to_index.size();
      if (npeff_fin < npeff) {
        cout << "Step " << it << ": Reducing grains from " << npeff << " to "
             << npeff_fin << endl;
        vector<double> Phi_shrink(pad_x * pad_y * pad_z * npeff_fin, 0.0);

        for (int ip = 0; ip < npeff_fin; ++ip) {
          int old_ip = index_to_index[ip];
          for (int iz = 0; iz < pad_z; ++iz) {
            for (int iy = 0; iy < pad_y; ++iy) {
              for (int ix = 0; ix < pad_x; ++ix) {
                Phi_shrink[idx4(ix, iy, iz, ip, npeff_fin)] =
                    Phi[idx4(ix, iy, iz, old_ip, npeff)];
              }
            }
          }
        }
        Phi = Phi_shrink;
        Phi_new = Phi;
        npeff = npeff_fin;
      }
    }
    if (it % 500 == 0) {
        cout << "Step " << it << " completed." << endl;
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

  for (int iz = 0; iz < pad_z; ++iz) {
    for (int iy = 0; iy < pad_y; ++iy) {
      for (int ix = 0; ix < pad_x; ++ix) {
        double local_sum = 0.0;
        for (int ip = 0; ip < npeff; ++ip)
          local_sum += Phi[idx4(ix, iy, iz, ip, npeff)];
        if (local_sum > 0) {
          for (int ip = 0; ip < npeff; ++ip)
            Phi[idx4(ix, iy, iz, ip, npeff)] /= local_sum;
        }
      }
    }
  }

  // Recalculate smsq (Cross terms only for grain boundaries)
  fill(smsq.begin(), smsq.end(), 0.0);
  for (int iz = 1; iz <= nz; ++iz) {
    for (int iy = 1; iy <= ny; ++iy) {
      for (int ix = 1; ix <= nx; ++ix) {
        double sum_cross = 0.0;
        for (int ip = 0; ip < npeff; ++ip) {
          for (int jp = 0; jp < npeff; ++jp) {
            if (ip != jp) {
              sum_cross += Phi[idx4(ix, iy, iz, ip, npeff)] *
                           Phi[idx4(ix, iy, iz, jp, npeff)];
            }
          }
        }
        smsq[idx3(ix, iy, iz)] = sum_cross;
      }
    }
  }

  // --- FILE EXPORT ---
  // Extract interior grid for outputs (ignoring pad layers, mapping to exact
  // dimensions make_props expects)

  cout << "writing into " << innerx << "x" << innery << "x" << innerz << endl;

  // 1. grain_vis.dat
  vector<double> smsq_out(innerx * innery * innerz);
  for (int iz = 0; iz <= innerz-1; ++iz) {
    for (int iy = 0; iy <= innery-1; ++iy) {
      for (int ix = 0; ix <= innerx-1; ++ix) {
        smsq_out[idx3_small(ix,iy,iz)] =
            smsq[idx3(ix+extrax, iy+extrax, iz+extrax)];
      }
    }
  }

  ofstream outVis("grain_vis.dat", ios::binary);
  if (outVis.is_open()) {
    outVis.write(reinterpret_cast<char *>(smsq_out.data()),
                 smsq_out.size() * sizeof(double));
    outVis.close();
  }

  // 2. grains.dat (This must perfectly match the 4D interior space for your
  // other C++ script)
  vector<double> Phi_out(innerx * innery * innerz * npeff);
  for (int iz = 0; iz <= innerz-1; ++iz) {
    for (int iy = 0; iy <= innery-1; ++iy) {
      for (int ix = 0; ix <= innerx-1; ++ix) {
        for (int ip = 0; ip < npeff; ++ip) {
          Phi_out[idx4_small(ix,iy,iz,ip, npeff)] = Phi[idx4(ix+extrax, iy+extrax, iz+extrax, ip, npeff)];
        }
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
