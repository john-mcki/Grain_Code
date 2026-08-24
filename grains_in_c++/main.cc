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
using Mat33 = array<array<double, 3>, 3>;
using Tensor4 = array<array<array<array<double, 3>, 3>, 3>, 3>;
const double PI = acos(-1.0);

// --- TENSOR MODULE EQUIVALENT ---
namespace TensorMod {
Tensor4 Cel_tnsr;

// Helper for 3x3 matrix multiplication
Mat33 matmul(const Mat33 &A, const Mat33 &B) {
  Mat33 C = {{{0}}};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 3; ++k)
        C[i][j] += A[i][k] * B[k][j];
  return C;
}

// Helper for 3x3 matrix transpose
Mat33 transpose(const Mat33 &A) {
  Mat33 C;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      C[i][j] = A[j][i];
  return C;
}

Mat33 RotX(double ang) {
  Mat33 R = {{{0}}};
  R[0][0] = 1.0;
  R[1][1] = cos(ang);
  R[1][2] = -sin(ang);
  R[2][1] = sin(ang);
  R[2][2] = cos(ang);
  return R;
}

Mat33 RotZ(double ang) {
  Mat33 R = {{{0}}};
  R[0][0] = cos(ang);
  R[0][1] = -sin(ang);
  R[1][0] = sin(ang);
  R[1][1] = cos(ang);
  R[2][2] = 1.0;
  return R;
}

void load_Cel_tnsr(const vector<double> &Cel_vec_in) {
  array<array<double, 6>, 6> Cel_mtrx = {{{0}}};

  // 0-indexed mapping of Voigt notation
  Cel_mtrx[0][0] = Cel_vec_in[0];
  Cel_mtrx[1][1] = Cel_vec_in[1];
  Cel_mtrx[2][2] = Cel_vec_in[2];
  Cel_mtrx[3][3] = Cel_vec_in[3];
  Cel_mtrx[4][4] = Cel_vec_in[4];
  Cel_mtrx[5][5] = Cel_vec_in[5];
  Cel_mtrx[0][1] = Cel_vec_in[6];
  Cel_mtrx[0][2] = Cel_vec_in[7];
  Cel_mtrx[0][3] = Cel_vec_in[8];
  Cel_mtrx[0][4] = Cel_vec_in[9];
  Cel_mtrx[0][5] = Cel_vec_in[10];
  Cel_mtrx[1][2] = Cel_vec_in[11];
  Cel_mtrx[1][3] = Cel_vec_in[12];
  Cel_mtrx[1][4] = Cel_vec_in[13];
  Cel_mtrx[1][5] = Cel_vec_in[14];
  Cel_mtrx[2][3] = Cel_vec_in[15];
  Cel_mtrx[2][4] = Cel_vec_in[16];
  Cel_mtrx[2][5] = Cel_vec_in[17];
  Cel_mtrx[3][4] = Cel_vec_in[18];
  Cel_mtrx[3][5] = Cel_vec_in[19];
  Cel_mtrx[4][5] = Cel_vec_in[20];

  // Symmetrize
  for (int i = 1; i < 6; ++i)
    for (int j = 0; j < i; ++j)
      Cel_mtrx[i][j] = Cel_mtrx[j][i];

  // Convert to 3x3x3x3
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        for (int l = 0; l < 3; ++l) {
          int a = (i == j) ? i : ((i + j == 1) ? 5 : ((i + j == 2) ? 4 : 3));
          int b = (k == l) ? k : ((k + l == 1) ? 5 : ((k + l == 2) ? 4 : 3));
          Cel_tnsr[i][j][k][l] = Cel_mtrx[a][b];
        }
      }
    }
  }
}

Mat33 load_Bel_tnsr(const vector<double> &Bel_vec_in) {
  Mat33 Bel_tnsr = {{{0}}};
  Bel_tnsr[0][0] = Bel_vec_in[0];
  Bel_tnsr[1][1] = Bel_vec_in[1];
  Bel_tnsr[2][2] = Bel_vec_in[2];
  Bel_tnsr[1][2] = Bel_vec_in[3];
  Bel_tnsr[0][2] = Bel_vec_in[4];
  Bel_tnsr[0][1] = Bel_vec_in[5];

  for (int i = 1; i < 3; ++i)
    for (int j = 0; j < i; ++j)
      Bel_tnsr[i][j] = Bel_tnsr[j][i];
  return Bel_tnsr;
}

vector<double> rot_Cel_tnsr(double ang1, double ang2, double ang3) {
  Mat33 Q = transpose(matmul(RotZ(ang1), matmul(RotX(ang2), RotZ(ang3))));
  Tensor4 Cel_tnsr_rot = {0};

  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 3; ++k)
        for (int l = 0; l < 3; ++l)
          for (int m = 0; m < 3; ++m)
            for (int n = 0; n < 3; ++n)
              for (int o = 0; o < 3; ++o)
                for (int p = 0; p < 3; ++p)
                  Cel_tnsr_rot[i][j][k][l] += Q[i][m] * Q[j][n] * Q[k][o] *
                                              Q[l][p] * Cel_tnsr[m][n][o][p];

  vector<double> Cel_vec_rot(21, 0.0);
  // Explicitly map back to Voigt vector (0-indexed)
  Cel_vec_rot[0] = Cel_tnsr_rot[0][0][0][0];
  Cel_vec_rot[1] = Cel_tnsr_rot[1][1][1][1];
  Cel_vec_rot[2] = Cel_tnsr_rot[2][2][2][2];
  Cel_vec_rot[3] = Cel_tnsr_rot[1][2][1][2];
  Cel_vec_rot[4] = Cel_tnsr_rot[0][2][0][2];
  Cel_vec_rot[5] = Cel_tnsr_rot[0][1][0][1];
  Cel_vec_rot[6] = Cel_tnsr_rot[0][0][1][1];
  Cel_vec_rot[7] = Cel_tnsr_rot[0][0][2][2];
  Cel_vec_rot[8] = Cel_tnsr_rot[0][0][1][2];
  Cel_vec_rot[9] = Cel_tnsr_rot[0][0][0][2];
  Cel_vec_rot[10] = Cel_tnsr_rot[0][0][0][1];
  Cel_vec_rot[11] = Cel_tnsr_rot[1][1][2][2];
  Cel_vec_rot[12] = Cel_tnsr_rot[1][1][1][2];
  Cel_vec_rot[13] = Cel_tnsr_rot[1][1][0][2];
  Cel_vec_rot[14] = Cel_tnsr_rot[1][1][0][1];
  Cel_vec_rot[15] = Cel_tnsr_rot[2][2][1][2];
  Cel_vec_rot[16] = Cel_tnsr_rot[2][2][0][2];
  Cel_vec_rot[17] = Cel_tnsr_rot[2][2][0][1];
  Cel_vec_rot[18] = Cel_tnsr_rot[1][2][0][2];
  Cel_vec_rot[19] = Cel_tnsr_rot[1][2][0][1];
  Cel_vec_rot[20] = Cel_tnsr_rot[0][2][0][1];

  // Precision cleanup
  double max_val = 0;
  for (double v : Cel_vec_rot)
    max_val = max(max_val, abs(v));
  for (double &v : Cel_vec_rot)
    if (abs(v) <= 1e10 * 2.22045e-16 * max_val)
      v = 0.0;

  return Cel_vec_rot;
}

vector<double> rot_Bel_tnsr(double ang1, double ang2, double ang3,
                            const Mat33 &Bel_tnsr) {
  Mat33 Q = transpose(matmul(RotZ(ang1), matmul(RotX(ang2), RotZ(ang3))));
  Mat33 Bel_tnsr_rot = {{{0}}};

  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      for (int m = 0; m < 3; ++m)
        for (int n = 0; n < 3; ++n)
          Bel_tnsr_rot[i][j] += Q[i][m] * Q[j][n] * Bel_tnsr[m][n];

  vector<double> Bel_vec_rot(6, 0.0);
  Bel_vec_rot[0] = Bel_tnsr_rot[0][0];
  Bel_vec_rot[1] = Bel_tnsr_rot[1][1];
  Bel_vec_rot[2] = Bel_tnsr_rot[2][2];
  Bel_vec_rot[3] = Bel_tnsr_rot[1][2];
  Bel_vec_rot[4] = Bel_tnsr_rot[0][2];
  Bel_vec_rot[5] = Bel_tnsr_rot[0][1];

  return Bel_vec_rot;
}
} // namespace TensorMod

// --- REDUCTION FUNCTIONS ---
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

// --- MAIN PROGRAM ---
int main() {
  // Dimensions
  const int nx = 256, ny = 256, nz = 1, np = 15;
  const int pad_x = nx + 2, pad_y = ny + 2, pad_z = nz + 2;

  // Initialize RNG
  mt19937 gen(2000);
  uniform_real_distribution<> dis(0.0, 1.0);

  vector<double> euler1(np), euler2(np), euler3(np);

  for (int i = 0; i < np; ++i) {
    vector<double> v1(3), v2(3), v3(3);
    double n1, n2;

    do {
      v1 = {dis(gen), dis(gen), dis(gen)};
      n1 = sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);
    } while (n1 == 0);
    for (int k = 0; k < 3; ++k)
      v1[k] /= n1;

    while (true) {
      do {
        v2 = {dis(gen), dis(gen), dis(gen)};
        n2 = sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]);
      } while (n2 == 0);
      for (int k = 0; k < 3; ++k)
        v2[k] /= n2;

      double dot = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
      if (abs(dot) < 0.99) {
        for (int k = 0; k < 3; ++k)
          v2[k] = v2[k] - dot * v1[k];
        n2 = sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]);
        for (int k = 0; k < 3; ++k)
          v2[k] /= n2;
        break;
      }
    }

    v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
    v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
    v3[2] = v1[0] * v2[1] - v1[1] * v2[0];

    euler1[i] = atan2(v3[1], v3[0]);
    euler2[i] = atan2(v3[0] * cos(euler1[i]) + v3[1] * sin(euler1[i]), v3[2]);
    euler3[i] = atan2(v2[2], -v1[2]);
  }

  vector<double> Cel_ref = {260, 260, 200, 45, 45, 94, 84, 46, -16, 0,  0,
                            46,  16,  0,   0,  0,  0,  0,  0,  0,   -14};
  vector<double> Bel_ref = {0.028, 0.028, -0.030, 0, 0, 0};
  vector<double> D_ref = {1e-3, 1e-3, 1e-4, 0, 0, 0};

  TensorMod::load_Cel_tnsr(Cel_ref);
  Mat33 Bel_tnsr = TensorMod::load_Bel_tnsr(Bel_ref);
  Mat33 D_tnsr = TensorMod::load_Bel_tnsr(D_ref);

  vector<vector<double>> CelRot(np, vector<double>(21));
  vector<vector<double>> BelRot(np, vector<double>(6));
  vector<vector<double>> DRot(np, vector<double>(6));

  for (int i = 0; i < np; ++i) {
    CelRot[i] = TensorMod::rot_Cel_tnsr(euler1[i], euler2[i], euler3[i]);
    BelRot[i] =
        TensorMod::rot_Bel_tnsr(euler1[i], euler2[i], euler3[i], Bel_tnsr);
    DRot[i] = TensorMod::rot_Bel_tnsr(euler1[i], euler2[i], euler3[i], D_tnsr);
  }

  // Allocate memory for spatial arrays
  // Layout: [x][y][z][p], flattened for efficiency
  auto idx3D = [pad_y, pad_z](int x, int y, int z) {
    return (x * pad_y * pad_z) + (y * pad_z) + z;
  };
  auto idx4D = [pad_y, pad_z, np](int x, int y, int z, int p) {
    return ((x * pad_y * pad_z) + (y * pad_z) + z) * np + p;
  };

  vector<double> PhiSolid(pad_x * pad_y * pad_z * np, 0.0);

  // Read grains.dat (assuming C++ binary formatting matches Fortran)
  cout << "reading in structure" << endl;

  // 1. Read into an unpadded temporary array
  vector<double> raw_grains(nx * ny * nz * np, 0.0);
  ifstream grainsFile("grains.dat", ios::binary);
  if (grainsFile.is_open()) {
    grainsFile.read(reinterpret_cast<char *>(raw_grains.data()),
                    raw_grains.size() * sizeof(double));
    grainsFile.close();

    // 2. Map the unpadded data into the interior of the padded PhiSolid array
    for (int i = 0; i < nx; ++i) {
      for (int j = 0; j < ny; ++j) {
        for (int k = 0; k < nz; ++k) {
          for (int p = 0; p < np; ++p) {
            int raw_idx = ((i * ny + j) * nz + k) * np + p;
            PhiSolid[idx4D(i + 1, j + 1, k + 1, p)] = raw_grains[raw_idx];
          }
        }
      }
    }

    // 3. Apply Periodic Boundary Conditions to the PhiSolid buffer layers
    // This ensures elastic properties wrap correctly on the edges
    for (int p = 0; p < np; ++p) {
      for (int k = 1; k <= nz; ++k) {
        for (int j = 1; j <= ny; ++j) {
          PhiSolid[idx4D(0, j, k, p)] = PhiSolid[idx4D(nx, j, k, p)];
          PhiSolid[idx4D(nx + 1, j, k, p)] = PhiSolid[idx4D(1, j, k, p)];
        }
      }
      for (int i = 0; i <= nx + 1; ++i) {
        for (int k = 1; k <= nz; ++k) {
          PhiSolid[idx4D(i, 0, k, p)] = PhiSolid[idx4D(i, ny, k, p)];
          PhiSolid[idx4D(i, ny + 1, k, p)] = PhiSolid[idx4D(i, 1, k, p)];
        }
      }
      for (int i = 0; i <= nx + 1; ++i) {
        for (int j = 0; j <= ny + 1; ++j) {
          PhiSolid[idx4D(i, j, 0, p)] = PhiSolid[idx4D(i, j, nz, p)];
          PhiSolid[idx4D(i, j, nz + 1, p)] = PhiSolid[idx4D(i, j, 1, p)];
        }
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
    double current_phi = PhiSolid[idx4D(1, 1, 1, p)];
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
        PsiOut_flat.push_back(1.0 - PhiSolid[idx4D(i, j, 1, outgrain)]);
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
  for (int k = 0; k < pad_z; ++k) {
    for (int j = 0; j < pad_y; ++j) {
      for (int i = 0; i < pad_x; ++i) {
        PhiSolid[idx4D(i, j, k, outgrain)] = 0.0;
      }
    }
  }

  cout << "renormalizing order parameters" << endl;
  for (int k = 0; k < pad_z; ++k) {
    for (int j = 0; j < pad_y; ++j) {
      for (int i = 0; i < pad_x; ++i) {
        double phisum = 0.0;
        for (int l = 0; l < np; ++l) {
          phisum += PhiSolid[idx4D(i, j, k, l)];
        }

        phisum = max(phisum, 1e-6);

        for (int l = 0; l < np; ++l) {
          PhiSolid[idx4D(i, j, k, l)] /= phisum;
        }
      }
    }
  }

  cout << "starting elastic constants calculation" << endl;

  // Indexing lambdas for the component-major 3D arrays
  auto tIdx21 = [pad_x, pad_y, pad_z](int c, int x, int y, int z) {
    return c * (pad_x * pad_y * pad_z) + (x * pad_y * pad_z) + (y * pad_z) + z;
  };
  auto tIdx6 = [pad_x, pad_y, pad_z](int c, int x, int y, int z) {
    return c * (pad_x * pad_y * pad_z) + (x * pad_y * pad_z) + (y * pad_z) + z;
  };

  // Allocate all three spatial fields
  vector<double> Cel(21 * pad_x * pad_y * pad_z, 0.0);
  vector<double> Bel(6 * pad_x * pad_y * pad_z, 0.0);
  vector<double> D_fld(6 * pad_x * pad_y * pad_z, 0.0);

  for (int l = 0; l < np; ++l) {
    for (int k = 0; k < pad_z; ++k) {
      for (int j = 0; j < pad_y; ++j) {
        for (int i = 0; i < pad_x; ++i) {
          double phi = PhiSolid[idx4D(i, j, k, l)];
          if (phi > 0.0) {
            for (int c = 0; c < 21; ++c)
              Cel[tIdx21(c, i, j, k)] += phi * CelRot[l][c];
            for (int c = 0; c < 6; ++c)
              Bel[tIdx6(c, i, j, k)] += phi * BelRot[l][c];
            for (int c = 0; c < 6; ++c)
              D_fld[tIdx6(c, i, j, k)] += phi * DRot[l][c];
          }
        }
      }
    }
  }

  // Apply Periodic Boundary Conditions to Cel (Matching Fortran)
  for (int c = 0; c < 21; ++c) {
    for (int j = 1; j <= ny + 1; ++j) {
      for (int k = 1; k <= nz + 1; ++k) {
        Cel[tIdx21(c, nx + 1, j, k)] = Cel[tIdx21(c, 1, j, k)];
      }
    }
    for (int i = 1; i <= nx + 1; ++i) {
      for (int k = 1; k <= nz + 1; ++k) {
        Cel[tIdx21(c, i, ny + 1, k)] = Cel[tIdx21(c, i, 1, k)];
      }
    }
    for (int i = 1; i <= nx + 1; ++i) {
      for (int j = 1; j <= ny + 1; ++j) {
        Cel[tIdx21(c, i, j, nz + 1)] = Cel[tIdx21(c, i, j, 1)];
      }
    }
  }

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

  // --- FILE EXPORT ---
  // Helper to extract specific tensor components into a flat array using
  // Fortran column-major order
  auto write_binary_subset = [&](const string &filename,
                                 const vector<double> &arr, int num_comps,
                                 int c_start, int c_end) {
    vector<double> out_data;
    out_data.reserve((c_end - c_start + 1) * (nx + 1) * (ny + 1));

    for (int iy = 1; iy <= ny + 1; ++iy) {
      for (int ix = 1; ix <= nx + 1; ++ix) {
        for (int c = c_start; c <= c_end; ++c) {
          int idx = (ix * pad_y * pad_z + iy * pad_z + 1) * num_comps + c;
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
  write_binary_subset("Cv2D_1.dat", Cel2D, 6, 0, 1);
  write_binary_subset("Cv2D_2.dat", Cel2D, 6, 2, 3);
  write_binary_subset("Cv2D_3.dat", Cel2D, 6, 4, 5);

  write_binary_subset("eigv2D_1.dat", Bel2D, 3, 0, 1);
  write_binary_subset("eigv2D_2.dat", Bel2D, 3, 2, 2);

  write_binary_subset("Dv2D_1.dat", D2D, 3, 0, 1);
  write_binary_subset("Dv2D_2.dat", D2D, 3, 2, 2);

  return 0;
}