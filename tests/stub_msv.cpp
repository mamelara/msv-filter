/*******************************************************************************
 * File: tests/stub_msv.cpp
 * Description: Stub implementation of compute_msv for compilation.
 * This provides a working MSV algorithm implementation for testing.
 * User should replace this with their optimized implementation.
 ******************************************************************************/

#include "hmmer_types.hpp"
#include "profile.hpp"
#include "dp_matrix.hpp"
#include <algorithm>
#include <cmath>

// MSV (Multiple Segment Viterbi) implementation
// Finds the maximum scoring segment alignment between a sequence and an HMM profile
//
// Algorithm:
// 1. Initialize DP[0][k] = 0 for all k (no score before sequence starts)
// 2. For each sequence position i from 1 to L:
//    For each model position k from 1 to M:
//      DP[i][k] = max(
//          DP[i-1][k-1] + score(i,k),  // extend from previous match
//          score(i,k)                   // start new segment here
//      )
// 3. Return max(DP[i][k]) over all i,k
//
// This is a simplified ungapped MSV - the real HMMER MSV is more complex
float compute_msv(const DigitalResidue* digital_sequence, int sequence_length,
                  const HMMProfile& profile, DPMatrix& dp_matrix, float expected_hit_count) {
    
    // Handle edge cases
    if (sequence_length <= 0 || profile.model_length <= 0) {
        return 0.0f;
    }
    
    const int M = profile.model_length;
    const int L = sequence_length;
    
    // Initialize first row of DP matrix (i=0, before sequence starts)
    for (int k = 0; k <= M; k++) {
        dp_matrix.match(0, k) = 0.0f;  // No score before sequence
    }
    
    float max_score = 0.0f;
    
    // Fill DP matrix
    for (int i = 1; i <= L; i++) {
        DigitalResidue residue = digital_sequence[i];
        
        // Skip invalid residues
        if (residue >= 20) {
            for (int k = 1; k <= M; k++) {
                dp_matrix.match(i, k) = 0.0f;
            }
            continue;
        }
        
        for (int k = 1; k <= M; k++) {
            float match_score = profile.match_score(k, residue);
            
            // MSV recurrence:
            // Option 1: Start a new segment at this position
            float start_new = match_score;
            
            // Option 2: Extend the previous segment
            // Note: For strict ungapped MSV, we only look at (i-1, k-1)
            float extend_prev = dp_matrix.match(i - 1, k - 1) + match_score;
            
            // Take the maximum
            float dp_value = std::max(start_new, extend_prev);
            
            // Ensure non-negative (MSV finds positive-scoring segments)
            dp_value = std::max(0.0f, dp_value);
            
            dp_matrix.match(i, k) = dp_value;
            
            // Track global maximum
            if (dp_value > max_score) {
                max_score = dp_value;
            }
        }
    }
    
    // Return the maximum scoring segment found
    // If all scores are negative, max_score will remain 0 (empty segment)
    return max_score;
}
