#include <btas/DENSE/Dreindex.h>

template<>
void btas::Dreindex<1>(const double* x, double* y, const btas::IVector<1>& xstr, const btas::IVector<1>& yshape)
{
  int i = 0;
  for(int i0 = 0; i0 < yshape[0]; ++i0) {
    y[i++] = x[i0];
  }
}

template<>
void btas::Dreindex<2>(const double* x, double* y, const btas::IVector<2>& xstr, const btas::IVector<2>& yshape)
{
  int i = 0;
  for(int i0 = 0; i0 < yshape[0]; ++i0) {
    int j0 = i0 * xstr[0];
    for(int i1 = 0; i1 < yshape[1]; ++i1) {
      y[i++] = x[j0 + i1 * xstr[1]];
    }
  }
}

template<>
void btas::Dreindex<3>(const double* x, double* y, const btas::IVector<3>& xstr, const btas::IVector<3>& yshape)
{
  int i = 0;
  for(int i0 = 0; i0 < yshape[0]; ++i0) {
    int j0 = i0 * xstr[0];
    for(int i1 = 0; i1 < yshape[1]; ++i1) {
      int j1 = j0 + i1 * xstr[1];
      for(int i2 = 0; i2 < yshape[2]; ++i2) {
        y[i++] = x[j1 + i2 * xstr[2]];
      }
    }
  }
}

template<>
void btas::Dreindex<4>(const double* x, double* y, const btas::IVector<4>& xstr, const btas::IVector<4>& yshape)
{
  int i = 0;
  for(int i0 = 0; i0 < yshape[0]; ++i0) {
    int j0 = i0 * xstr[0];
    for(int i1 = 0; i1 < yshape[1]; ++i1) {
      int j1 = j0 + i1 * xstr[1];
      for(int i2 = 0; i2 < yshape[2]; ++i2) {
        int j2 = j1 + i2 * xstr[2];
        for(int i3 = 0; i3 < yshape[3]; ++i3) {
          y[i++] = x[j2 + i3 * xstr[3]];
        }
      }
    }
  }
}

template<>
void btas::Dreindex<5>(const double* x, double* y, const btas::IVector<5>& xstr, const btas::IVector<5>& yshape)
{
  int i = 0;
  for(int i0 = 0; i0 < yshape[0]; ++i0) {
    int j0 = i0 * xstr[0];
    for(int i1 = 0; i1 < yshape[1]; ++i1) {
      int j1 = j0 + i1 * xstr[1];
      for(int i2 = 0; i2 < yshape[2]; ++i2) {
        int j2 = j1 + i2 * xstr[2];
        for(int i3 = 0; i3 < yshape[3]; ++i3) {
          int j3 = j2 + i3 * xstr[3];
          for(int i4 = 0; i4 < yshape[4]; ++i4) {
            y[i++] = x[j3 + i4 * xstr[4]];
          }
        }
      }
    }
  }
}

template<>
void btas::Dreindex<6>(const double* x, double* y, const btas::IVector<6>& xstr, const btas::IVector<6>& yshape)
{
  int i = 0;
  for(int i0 = 0; i0 < yshape[0]; ++i0) {
    int j0 = i0 * xstr[0];
    for(int i1 = 0; i1 < yshape[1]; ++i1) {
      int j1 = j0 + i1 * xstr[1];
      for(int i2 = 0; i2 < yshape[2]; ++i2) {
        int j2 = j1 + i2 * xstr[2];
        for(int i3 = 0; i3 < yshape[3]; ++i3) {
          int j3 = j2 + i3 * xstr[3];
          for(int i4 = 0; i4 < yshape[4]; ++i4) {
            int j4 = j3 + i4 * xstr[4];
            for(int i5 = 0; i5 < yshape[5]; ++i5) {
              y[i++] = x[j4 + i5 * xstr[5]];
            }
          }
        }
      }
    }
  }
}

template<>
void btas::Dreindex<7>(const double* x, double* y, const btas::IVector<7>& xstr, const btas::IVector<7>& yshape)
{
  int i = 0;
  for(int i0 = 0; i0 < yshape[0]; ++i0) {
    int j0 = i0 * xstr[0];
    for(int i1 = 0; i1 < yshape[1]; ++i1) {
      int j1 = j0 + i1 * xstr[1];
      for(int i2 = 0; i2 < yshape[2]; ++i2) {
        int j2 = j1 + i2 * xstr[2];
        for(int i3 = 0; i3 < yshape[3]; ++i3) {
          int j3 = j2 + i3 * xstr[3];
          for(int i4 = 0; i4 < yshape[4]; ++i4) {
            int j4 = j3 + i4 * xstr[4];
            for(int i5 = 0; i5 < yshape[5]; ++i5) {
              int j5 = j4 + i5 * xstr[5];
              for(int i6 = 0; i6 < yshape[6]; ++i6) {
                y[i++] = x[j5 + i6 * xstr[6]];
              }
            }
          }
        }
      }
    }
  }
}

template<>
void btas::Dreindex<8>(const double* x, double* y, const btas::IVector<8>& xstr, const btas::IVector<8>& yshape)
{
  int i = 0;
  for(int i0 = 0; i0 < yshape[0]; ++i0) {
    int j0 = i0 * xstr[0];
    for(int i1 = 0; i1 < yshape[1]; ++i1) {
      int j1 = j0 + i1 * xstr[1];
      for(int i2 = 0; i2 < yshape[2]; ++i2) {
        int j2 = j1 + i2 * xstr[2];
        for(int i3 = 0; i3 < yshape[3]; ++i3) {
          int j3 = j2 + i3 * xstr[3];
          for(int i4 = 0; i4 < yshape[4]; ++i4) {
            int j4 = j3 + i4 * xstr[4];
            for(int i5 = 0; i5 < yshape[5]; ++i5) {
              int j5 = j4 + i5 * xstr[5];
              for(int i6 = 0; i6 < yshape[6]; ++i6) {
                int j6 = j5 + i6 * xstr[6];
                for(int i7 = 0; i7 < yshape[7]; ++i7) {
                  y[i++] = x[j6 + i7 * xstr[7]];
                }
              }
            }
          }
        }
      }
    }
  }
}

template<>
void btas::Dreindex<9>(const double* x, double* y, const btas::IVector<9>& xstr, const btas::IVector<9>& yshape)
{
  int i = 0;
  for(int i0 = 0; i0 < yshape[0]; ++i0) {
    int j0 = i0 * xstr[0];
    for(int i1 = 0; i1 < yshape[1]; ++i1) {
      int j1 = j0 + i1 * xstr[1];
      for(int i2 = 0; i2 < yshape[2]; ++i2) {
        int j2 = j1 + i2 * xstr[2];
        for(int i3 = 0; i3 < yshape[3]; ++i3) {
          int j3 = j2 + i3 * xstr[3];
          for(int i4 = 0; i4 < yshape[4]; ++i4) {
            int j4 = j3 + i4 * xstr[4];
            for(int i5 = 0; i5 < yshape[5]; ++i5) {
              int j5 = j4 + i5 * xstr[5];
              for(int i6 = 0; i6 < yshape[6]; ++i6) {
                int j6 = j5 + i6 * xstr[6];
                for(int i7 = 0; i7 < yshape[7]; ++i7) {
                  int j7 = j6 + i7 * xstr[7];
                  for(int i8 = 0; i8 < yshape[8]; ++i8) {
                    y[i++] = x[j7 + i8 * xstr[8]];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

template<>
void btas::Dreindex<10>(const double* x, double* y, const btas::IVector<10>& xstr, const btas::IVector<10>& yshape)
{
  int i = 0;
  for(int i0 = 0; i0 < yshape[0]; ++i0) {
    int j0 = i0 * xstr[0];
    for(int i1 = 0; i1 < yshape[1]; ++i1) {
      int j1 = j0 + i1 * xstr[1];
      for(int i2 = 0; i2 < yshape[2]; ++i2) {
        int j2 = j1 + i2 * xstr[2];
        for(int i3 = 0; i3 < yshape[3]; ++i3) {
          int j3 = j2 + i3 * xstr[3];
          for(int i4 = 0; i4 < yshape[4]; ++i4) {
            int j4 = j3 + i4 * xstr[4];
            for(int i5 = 0; i5 < yshape[5]; ++i5) {
              int j5 = j4 + i5 * xstr[5];
              for(int i6 = 0; i6 < yshape[6]; ++i6) {
                int j6 = j5 + i6 * xstr[6];
                for(int i7 = 0; i7 < yshape[7]; ++i7) {
                  int j7 = j6 + i7 * xstr[7];
                  for(int i8 = 0; i8 < yshape[8]; ++i8) {
                    int j8 = j7 + i8 * xstr[8];
                    for(int i9 = 0; i9 < yshape[9]; ++i9) {
                      y[i++] = x[j8 + i9 * xstr[9]];
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

template<>
void btas::Dreindex<11>(const double* x, double* y, const btas::IVector<11>& xstr, const btas::IVector<11>& yshape)
{
  int i = 0;
  for(int i0 = 0; i0 < yshape[0]; ++i0) {
    int j0 = i0 * xstr[0];
    for(int i1 = 0; i1 < yshape[1]; ++i1) {
      int j1 = j0 + i1 * xstr[1];
      for(int i2 = 0; i2 < yshape[2]; ++i2) {
        int j2 = j1 + i2 * xstr[2];
        for(int i3 = 0; i3 < yshape[3]; ++i3) {
          int j3 = j2 + i3 * xstr[3];
          for(int i4 = 0; i4 < yshape[4]; ++i4) {
            int j4 = j3 + i4 * xstr[4];
            for(int i5 = 0; i5 < yshape[5]; ++i5) {
              int j5 = j4 + i5 * xstr[5];
              for(int i6 = 0; i6 < yshape[6]; ++i6) {
                int j6 = j5 + i6 * xstr[6];
                for(int i7 = 0; i7 < yshape[7]; ++i7) {
                  int j7 = j6 + i7 * xstr[7];
                  for(int i8 = 0; i8 < yshape[8]; ++i8) {
                    int j8 = j7 + i8 * xstr[8];
                    for(int i9 = 0; i9 < yshape[9]; ++i9) {
                      int j9 = j8 + i9 * xstr[9];
                      for(int i10 = 0; i10 < yshape[10]; ++i10) {
                        y[i++] = x[j9 + i10 * xstr[10]];
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

template<>
void btas::Dreindex<12>(const double* x, double* y, const btas::IVector<12>& xstr, const btas::IVector<12>& yshape)
{
  int i = 0;
  for(int i0 = 0; i0 < yshape[0]; ++i0) {
    int j0 = i0 * xstr[0];
    for(int i1 = 0; i1 < yshape[1]; ++i1) {
      int j1 = j0 + i1 * xstr[1];
      for(int i2 = 0; i2 < yshape[2]; ++i2) {
        int j2 = j1 + i2 * xstr[2];
        for(int i3 = 0; i3 < yshape[3]; ++i3) {
          int j3 = j2 + i3 * xstr[3];
          for(int i4 = 0; i4 < yshape[4]; ++i4) {
            int j4 = j3 + i4 * xstr[4];
            for(int i5 = 0; i5 < yshape[5]; ++i5) {
              int j5 = j4 + i5 * xstr[5];
              for(int i6 = 0; i6 < yshape[6]; ++i6) {
                int j6 = j5 + i6 * xstr[6];
                for(int i7 = 0; i7 < yshape[7]; ++i7) {
                  int j7 = j6 + i7 * xstr[7];
                  for(int i8 = 0; i8 < yshape[8]; ++i8) {
                    int j8 = j7 + i8 * xstr[8];
                    for(int i9 = 0; i9 < yshape[9]; ++i9) {
                      int j9 = j8 + i9 * xstr[9];
                      for(int i10 = 0; i10 < yshape[10]; ++i10) {
                        int j10 = j9 + i10 * xstr[10];
                        for(int i11 = 0; i11 < yshape[11]; ++i11) {
                          y[i++] = x[j10 + i11 * xstr[11]];
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}


