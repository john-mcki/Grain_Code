#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

using namespace std;

// --- TYPE ALIASES & CONSTANTS ---
using Mat22 = array<array<double, 2>, 2>;
using Tensor4_2D = array<array<array<array<double, 2>, 2>, 2>, 2>;
const double PI = acos(-1.0);

// --- TENSOR MODULE EQUIVALENT ---
namespace TensorMod {
// Tensor4 Cel_tnsr;
Tensor4_2D Cel_tnsr_2D;

// Helper for 2x2 matrix multiplication
Mat22 matmul(const Mat22 &A, const Mat22 &B) {
  Mat22 C = {{{0}}};
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < 2; ++k)
        C[i][j] += A[i][k] * B[k][j];
  return C;
}

// Helper for 2x2 matrix transpose
Mat22 transpose(const Mat22 &A) {
  Mat22 C;
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
      C[i][j] = A[j][i];
  return C;
}

// 2x2 matrix rotation about the Z-axis
Mat22 RotZ_2D(double ang) {
  Mat22 R = {{{0}}};
  R[0][0] = cos(ang);
  R[0][1] = -sin(ang);
  R[1][0] = sin(ang);
  R[1][1] = cos(ang);
  return R;
}

void load_Cel_tnsr_2D(const vector<double> &Cel_vec) {
  array<array<double, 3>, 3> Cel_mtrx = {{{0}}};

  // 0-indexed mapping of Voigt notation for 2D
  Cel_mtrx[0][0] = Cel_vec[0]; // C11
  Cel_mtrx[1][1] = Cel_vec[1]; // C22
  Cel_mtrx[2][2] = Cel_vec[2]; // C66
  Cel_mtrx[0][1] = Cel_vec[3]; // C12
  Cel_mtrx[0][2] = Cel_vec[4]; // C16
  Cel_mtrx[1][2] = Cel_vec[5]; // C26

  // Symmetrize
  for (int i = 1; i < 3; ++i)
    for (int j = 0; j < i; ++j)
      Cel_mtrx[i][j] = Cel_mtrx[j][i];

  // Convert to 2x2x2x2
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      for (int k = 0; k < 2; ++k) {
        for (int l = 0; l < 2; ++l) {
          int a = (i == j) ? i : ((i + j == 1) ? 2 : 3);
          int b = (k == l) ? k : ((k + l == 1) ? 2 : 3);
          Cel_tnsr_2D[i][j][k][l] =
              Cel_mtrx[a][b]; // Full tensor needed for rotation
        }
      }
    }
  }
}

Mat22 load_Bel_tnsr_2D(const vector<double> &Bel_vec) {
  Mat22 Bel_tnsr = {{{0}}};
  Bel_tnsr[0][0] = Bel_vec[0];
  Bel_tnsr[1][1] = Bel_vec[1];
  Bel_tnsr[0][1] = Bel_vec[2];
  Bel_tnsr[1][0] = Bel_vec[2];
  return Bel_tnsr;
}

vector<double> rot_Cel_tnsr_2D(double ang) {
  Mat22 Q = transpose(RotZ_2D(ang));
  Tensor4_2D Cel_tnsr_rot = {0};

  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < 2; ++k)
        for (int l = 0; l < 2; ++l)
          for (int m = 0; m < 2; ++m)
            for (int n = 0; n < 2; ++n)
              for (int o = 0; o < 2; ++o)
                for (int p = 0; p < 2; ++p)
                  Cel_tnsr_rot[i][j][k][l] += Q[i][m] * Q[j][n] * Q[k][o] *
                                              Q[l][p] * Cel_tnsr_2D[m][n][o][p];

  vector<double> Cel_vec_rot(6, 0.0);
  Cel_vec_rot[0] = Cel_tnsr_rot[0][0][0][0]; // C11
  Cel_vec_rot[1] = Cel_tnsr_rot[1][1][1][1]; // C22
  Cel_vec_rot[2] = Cel_tnsr_rot[0][1][0][1]; // C66
  Cel_vec_rot[3] = Cel_tnsr_rot[0][0][1][1]; // C12
  Cel_vec_rot[4] = Cel_tnsr_rot[0][0][0][1]; // C16
  Cel_vec_rot[5] = Cel_tnsr_rot[1][1][0][1]; // C26

  // Precision cleanup
  double max_val = 0;
  for (double v : Cel_vec_rot)
    max_val = max(max_val, abs(v));
  for (double &v : Cel_vec_rot)
    if (abs(v) <= 1e10 * 2.22045e-16 * max_val)
      v = 0.0;

  return Cel_vec_rot;
}

vector<double> rot_Bel_tnsr_2D(double ang, const Mat22 &Bel_tnsr) {
  Mat22 Q = transpose(RotZ_2D(ang));
  Mat22 Bel_tnsr_rot = {{{0}}};

  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
      for (int m = 0; m < 2; ++m)
        for (int n = 0; n < 2; ++n)
          Bel_tnsr_rot[i][j] += Q[i][m] * Q[j][n] * Bel_tnsr[m][n];

  vector<double> Bel_vec_rot(3, 0.0);
  Bel_vec_rot[0] = Bel_tnsr_rot[0][0];
  Bel_vec_rot[1] = Bel_tnsr_rot[1][1];
  Bel_vec_rot[2] = Bel_tnsr_rot[0][1];

  return Bel_vec_rot;
}
} // namespace TensorMod

// --- REDUCTION FUNCTIONS --- //
void reduce2D_4thrank(const vector<double> &in, vector<double> &out,
                      int out_offset) {
  out[out_offset + 0] = in[0];  // C11
  out[out_offset + 1] = in[1];  // C22
  out[out_offset + 2] = in[5];  // C66
  out[out_offset + 3] = in[6];  // C12
  out[out_offset + 4] = in[10]; // C16
  out[out_offset + 5] = in[14]; // C26
}

void reduce2D_2ndrank(const vector<double> &in, vector<double> &out,
                      int out_offset) {
  out[out_offset + 0] = in[0];
  out[out_offset + 1] = in[1];
  out[out_offset + 2] = in[5];
}

// --- MAIN PROGRAM --- //
int main() {
  // Dimensions
  const int nx = 256, ny = 256, np = 29;
  const int pad_x = nx + 2, pad_y = ny + 2;

  // Initialize RNG
  mt19937 gen(2000);
  uniform_real_distribution<> dis(0.0, 1.0);

  // One theta angle for 2D rotation, not dependent on euler angles
  // Theta itself can be directly calculated from a randomly generated number
  // and multipled by pi, spectrum of possible theta values is thus (0, pi)
  // NOTE: follows normal distribution trend, possibly favored toward pi/2
  vector<double> theta(np);
  for (int i = 0; i < np; ++i) {
    theta[i] = dis(gen) * PI;
  }

  // Specific thetas for bicrystal case
  /*
  vector<double> theta(np);
  vector<double> bicrystal_rotations(np); // theta values
  bicrystal_rotations[0] = 0.0;
  bicrystal_rotations[1] = PI / 2.0;
  bicrystal_rotations[2] = 0.0;
  for (int i = 0; i < np; ++i) {
    theta[i] = bicrystal_rotations[i];
  }
  */

  // Check these later
  vector<double> Cel_ref_2D = {260, 260, 45, 84, -16, 0};
  vector<double> Bel_ref_2D = {0.028, 0.028, -0.030};
  vector<double> D_ref_2D = {1e-3, 1e-3, 1e-4};

  TensorMod::load_Cel_tnsr_2D(Cel_ref_2D);
  Mat22 Bel_tnsr_2D = TensorMod::load_Bel_tnsr_2D(Bel_ref_2D);
  Mat22 D_tnsr_2D = TensorMod::load_Bel_tnsr_2D(D_ref_2D);

  vector<vector<double>> CelRot_2D(np, vector<double>(6));
  vector<vector<double>> BelRot_2D(np, vector<double>(3));
  vector<vector<double>> DRot_2D(np, vector<double>(3));

  for (int i = 0; i < np; ++i) {
    CelRot_2D[i] = TensorMod::rot_Cel_tnsr_2D(theta[i]);
    BelRot_2D[i] = TensorMod::rot_Bel_tnsr_2D(theta[i], Bel_tnsr_2D);
    DRot_2D[i] = TensorMod::rot_Bel_tnsr_2D(theta[i], D_tnsr_2D);
  }

  auto idx2D = [pad_x, pad_y](int x, int y) { return (y * pad_x) + x; };

  auto idx3D = [pad_x, pad_y, np](int x, int y, int p) {
    return ((y * pad_x) + x) * np + p;
  };

  vector<double> PhiSolid(pad_x * pad_y * np, 0.0);

  // Read grains.dat (assuming C++ binary formatting matches Fortran)
  cout << "reading in structure" << endl;

  // 1. Read into an unpadded temporary array
  // vector<double> raw_grains(nx * ny * nz * np, 0.0);
  vector<double> raw_grains(nx * ny * np, 0.0);
  ifstream grainsFile("grains.dat", ios::binary);
  if (grainsFile.is_open()) {
    grainsFile.read(reinterpret_cast<char *>(raw_grains.data()),
                    raw_grains.size() * sizeof(double));
    grainsFile.close();

    // 2. Map the unpadded data into the interior of the padded PhiSolid array
    for (int i = 0; i < nx; ++i) {
      for (int j = 0; j < ny; ++j) {
        for (int p = 0; p < np; ++p) {
          int raw_idx = (i * ny + j) * np + p;
          PhiSolid[idx3D(i + 1, j + 1, p)] = raw_grains[raw_idx];
        }
      }
    }

    // 3. Apply Periodic Boundary Conditions to the PhiSolid buffer layers
    for (int p = 0; p < np; ++p) {
      for (int iy = 0; iy < pad_y; ++iy) {
        PhiSolid[idx3D(0, iy, p)] = PhiSolid[idx3D(nx, iy, p)];
        PhiSolid[idx3D(nx + 1, iy, p)] = PhiSolid[idx3D(1, iy, p)];
      }
      for (int ix = 0; ix < pad_x; ++ix) {
        PhiSolid[idx3D(ix, 0, p)] = PhiSolid[idx3D(ix, ny, p)];
        PhiSolid[idx3D(ix, ny + 1, p)] = PhiSolid[idx3D(ix, 1, p)];
      }
    }
  } else {
    cout << "Could not open grains.dat - Please ensure it exists." << endl;
  }

  // --- OUTER GRAIN LOGIC & PSI.DAT ---
  // 1. Find the outer grain ID (based on max volume fraction at coordinate
  // 1,1,1)
  int outgrain = 0;
  double max_phi = -1.0;
  for (int p = 0; p < np; ++p) {
    double current_phi = PhiSolid[idx3D(1, 1, p)];
    if (current_phi > max_phi) {
      max_phi = current_phi;
      outgrain = p;
    }
  }
  // Print 1-based ID to match Fortran output
  cout << "outer grainID: " << outgrain + 1 << endl;

  // 2. Generate and write PsiOut
  // Fortran array is PsiOut(x1:x2+1, y1:y2+1) which means bounds are nx+1 by
  // ny+1
  vector<double> PsiOut_flat;
  PsiOut_flat.reserve((nx + 1) * (ny + 1));

  // To match Fortran's column-major binary output, the X-axis (i) must be the
  // inner loop
  for (int j = 1; j <= ny + 1; ++j) {
    for (int i = 1; i <= nx + 1; ++i) {
      if (i <= nx && j <= ny) {
        PsiOut_flat.push_back(1.0 - PhiSolid[idx3D(i, j, outgrain)]);
      } else {
        PsiOut_flat.push_back(0.0);
      }
    }
  }

  ofstream psiFile("psi.dat", ios::binary);
  if (psiFile.is_open()) {
    psiFile.write(reinterpret_cast<char *>(PsiOut_flat.data()),
                  PsiOut_flat.size() * sizeof(double));
    psiFile.close();
  }

  // 3. Remove the outer grain and renormalize the remaining order parameters
  for (int j = 0; j < pad_y; ++j) {
    for (int i = 0; i < pad_x; ++i) {
      PhiSolid[idx3D(i, j, outgrain)] = 0.0;
    }
  }

  cout << "renormalizing order parameters" << endl;

  for (int j = 0; j < pad_y; ++j) {
    for (int i = 0; i < pad_x; ++i) {
      double phisum = 0.0;
      for (int l = 0; l < np; ++l) {
        phisum += PhiSolid[idx3D(i, j, l)];
      }

      phisum = max(phisum, 1e-6);

      for (int l = 0; l < np; ++l) {
        PhiSolid[idx3D(i, j, l)] /= phisum;
      }
    }
  }

  cout << "starting elastic constants calculation" << endl;

  auto tIdx6 = [pad_x, pad_y](int c, int x, int y) {
    return c * (pad_x * pad_y) + (x * pad_y) + y;
  };
  auto tIdx3 = [pad_x, pad_y](int c, int x, int y) {
    return c * (pad_x * pad_y) + (x * pad_y) + y;
  };

  // Allocate all two spatial fields
  vector<double> Cel(6 * pad_x * pad_y, 0.0);
  vector<double> Bel(3 * pad_x * pad_y, 0.0);
  vector<double> D_fld(3 * pad_x * pad_y, 0.0);

  for (int l = 0; l < np; ++l) {
    for (int j = 0; j < pad_y; ++j) {
      for (int i = 0; i < pad_x; ++i) {
        double phi = PhiSolid[idx3D(i, j, l)];
        if (phi > 0.0) {
          for (int c = 0; c < 6; ++c)
            Cel[tIdx6(c, i, j)] += phi * CelRot_2D[l][c];
          for (int c = 0; c < 3; ++c)
            Bel[tIdx3(c, i, j)] += phi * BelRot_2D[l][c];
          for (int c = 0; c < 3; ++c)
            D_fld[tIdx3(c, i, j)] += phi * DRot_2D[l][c];
        }
      }
    }
  }

  // Apply Periodic Boundary Conditions to Cel
  for (int c = 0; c < 6; ++c) {
    for (int j = 1; j <= ny + 1; ++j) {
      Cel[tIdx6(c, nx + 1, j)] = Cel[tIdx6(c, 1, j)];
    }
    for (int i = 1; i <= nx + 1; ++i) {
      Cel[tIdx6(c, i, ny + 1)] = Cel[tIdx6(c, i, 1)];
    }
  }

  // Fully unnecessary with my adjustments
  /*
  // Reduction to 2D
  vector<double> Cel2D(6 * pad_x * pad_y * pad_z, 0.0);
  vector<double> Bel2D(3 * pad_x * pad_y * pad_z, 0.0);
  vector<double> D2D(3 * pad_x * pad_y * pad_z, 0.0);

  for (int k = 0; k < pad_z; ++k) {
    for (int j = 0; j < pad_y; ++j) {
      for (int i = 0; i < pad_x; ++i) {
        vector<double> temp_Cel(21), temp_Bel(6), temp_D(6);
        for (int c = 0; c < 21; ++c)
          temp_Cel[c] = Cel[tIdx21(c, i, j, k)];
        for (int c = 0; c < 6; ++c)
          temp_Bel[c] = Bel[tIdx6(c, i, j, k)];
        for (int c = 0; c < 6; ++c)
          temp_D[c] = D_fld[tIdx6(c, i, j, k)];

        int off6 = (i * pad_y * pad_z + j * pad_z + k) * 6;
        int off3 = (i * pad_y * pad_z + j * pad_z + k) * 3;

        reduce2D_4thrank(temp_Cel, Cel2D, off6);
        reduce2D_2ndrank(temp_Bel, Bel2D, off3);
        reduce2D_2ndrank(temp_D, D2D, off3);
      }
    }
  }

  // Apply MAX limit to Cel2D components 1:3 (Indices 0, 1, 2)
  for (int iy = 1; iy <= ny + 1; ++iy) {
    for (int ix = 1; ix <= nx + 1; ++ix) {
      for (int c = 0; c < 3; ++c) {
        int idx = (ix * pad_y * pad_z + iy * pad_z + 1) * 6 + c;
        Cel2D[idx] = max(Cel2D[idx], 1.0);
      }
    }
  }
  */

  // --- FILE EXPORT ---
  // Helper to extract specific tensor components into a flat array using
  // Fortran column-major order
  /*
  auto write_binary_subset = [&](const string &filename,
                                 const vector<double> &arr, int num_comps,
                                 int c_start, int c_end) {
    vector<double> out_data;
    out_data.reserve((c_end - c_start + 1) * (nx + 1) * (ny + 1));

    for (int iy = 1; iy <= ny + 1; ++iy) {
      for (int ix = 1; ix <= nx + 1; ++ix) {
        for (int c = c_start; c <= c_end; ++c) {
          int idx = (ix * pad_y + iy + 1) * num_comps + c;
          out_data.push_back(arr[idx]);
        }
      }
    }

    ofstream out(filename, ios::binary);
    if (out.is_open()) {
      out.write(reinterpret_cast<char *>(out_data.data()),
                out_data.size() * sizeof(double));
    }
  };
  */

  // Helper to extract specific tensor components into a flat array using
  // Fortran column-major order
  auto write_binary_subset = [&](const string &filename,
                                 const vector<double> &arr, int c_start,
                                 int c_end) {
    vector<double> out_data;
    out_data.reserve((c_end - c_start + 1) * (nx + 1) * (ny + 1));

    for (int iy = 1; iy <= ny + 1; ++iy) {
      for (int ix = 1; ix <= nx + 1; ++ix) {
        for (int c = c_start; c <= c_end; ++c) {
          int idx = c * (pad_x * pad_y) + (ix * pad_y) + iy;
          out_data.push_back(arr[idx]);
        }
      }
    }

    ofstream out(filename, ios::binary);
    if (out.is_open()) {
      out.write(reinterpret_cast<char *>(out_data.data()),
                out_data.size() * sizeof(double));
    }
  };

  // Output identical data blocks to the 7 original Fortran files
  write_binary_subset("Cv2D_1.dat", Cel, 0, 1);
  write_binary_subset("Cv2D_2.dat", Cel, 2, 3);
  write_binary_subset("Cv2D_3.dat", Cel, 4, 5);

  write_binary_subset("eigv2D_1.dat", Bel, 0, 1);
  write_binary_subset("eigv2D_2.dat", Bel, 2, 2);

  write_binary_subset("Dv2D_1.dat", D_fld, 0, 1);
  write_binary_subset("Dv2D_2.dat", D_fld, 2, 2);

  return 0;
}