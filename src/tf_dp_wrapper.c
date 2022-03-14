#include "tf.h"
/**
 * @brief Wrapper called from R using .C
 * @param n                    number of observations
 * @param y                    response vector
 * @param lam                  the maximum lambda of the path
 * @param beta                 allocated space for the output
 * @return  void
 */
void tf_dp_wrapper(int* n, double* y, double* lam, double* beta)
{
    tf_dp(n[0], y, lam[0], beta);
}
