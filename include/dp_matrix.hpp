/*******************************************************************************
 * File: include/dp_matrix.hpp
 * Description: HMMER-compatible P7_GMX (generic DP matrix) structure for MSV
 * filter mocking. Replicates hmmer/src/hmmer.h
 ******************************************************************************/

#ifndef MSV_FILTER_DP_MATRIX_HPP
#define MSV_FILTER_DP_MATRIX_HPP

#include <vector>
#include <limits>
#include "hmmer_types.hpp"

/*******************************************************************************
 * P7_GMX Structure (from hmmer.h)
 * 
 * Generic Dynamic Programming Matrix
 * Used for storing DP values during MSV, Viterbi, and Forward algorithms
 ******************************************************************************/

class DPMatrix {
public:
    // --- Dimensions ---
    int model_length;  // Actual model dimension (model 1..model_length)
    int sequence_length;  // Actual sequence dimension (seq 1..sequence_length)
    
    // --- Allocation Info (for dynamic growth) ---
    int allocR;      // Allocated # of rows
    int validR;      // Valid # of rows
    int allocW;      // Allocated row width
    
    // --- DP Data ---
    
    // dp[i][k * p7G_NSCELLS + s] where:
    //   i = row index (0..sequence_length)
    //   k = model position (0..model_length)
    //   s = state (0=M, 1=I, 2=D)
    // MSV only uses s=0 (M state)
    std::vector<std::vector<float>> dp;
    
    // xmx[i * p7G_NXCELLS + s] where:
    //   i = row index (0..sequence_length)
    //   s = special state (0=E, 1=N, 2=J, 3=B, 4=C)
    std::vector<float> xmx;
    
    // --- Constructor ---
    DPMatrix(int max_model_length, int max_sequence_length)
        : model_length(max_model_length), sequence_length(max_sequence_length), 
          allocR(max_sequence_length + 1), validR(max_sequence_length + 1), allocW(max_model_length + 1)
    {
        // Allocate dp: (sequence_length+1) rows, each with (model_length+1) * 3 cells
        // Note: In HMMER, row 0 is for initialization (before sequence)
        dp.resize(sequence_length + 1);
        for (int i = 0; i <= sequence_length; i++) {
            dp[i].resize((model_length + 1) * p7G_NSCELLS, -eslINFINITY);
        }
        
        // Allocate xmx: (sequence_length+1) rows * 5 special states
        xmx.resize((sequence_length + 1) * p7G_NXCELLS, -eslINFINITY);
    }
    
    // --- Accessor Methods (replace HMMER macros) ---
    
    // MMX(i,k) = dp[(i)][(k) * p7G_NSCELLS + p7G_M]
    inline float& match(int i, int k) {
        return dp[i][(k * p7G_NSCELLS) + p7G_M];
    }
    
    inline float match(int i, int k) const {
        return dp[i][(k * p7G_NSCELLS) + p7G_M];
    }
    
    // IMX(i,k) = dp[(i)][(k) * p7G_NSCELLS + p7G_I]
    inline float& insert(int i, int k) {
        return dp[i][(k * p7G_NSCELLS) + p7G_I];
    }
    
    // DMX(i,k) = dp[(i)][(k) * p7G_NSCELLS + p7G_D]
    inline float& delete_state(int i, int k) {
        return dp[i][(k * p7G_NSCELLS) + p7G_D];
    }
    
    // XMX(i,s) = xmx[(i) * p7G_NXCELLS + (s)]
    inline float& special(int i, int s) {
        return xmx[(i * p7G_NXCELLS) + s];
    }
    
    inline float special(int i, int s) const {
        return xmx[(i * p7G_NXCELLS) + s];
    }
};

#endif // MSV_FILTER_DP_MATRIX_HPP
