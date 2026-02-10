/*******************************************************************************
 * File: include/hmmer_types.hpp
 * Description: HMMER-compatible type definitions and constants for MSV filter
 * mocking. These replicate the essential structures from hmmer/src/hmmer.h
 * and easel/easel.h
 ******************************************************************************/

#ifndef MSV_FILTER_HMMER_TYPES_HPP
#define MSV_FILTER_HMMER_TYPES_HPP

#include <cstdint>
#include <cfloat>
#include <cmath>
#include <vector>

/*******************************************************************************
 * 1. EASEL Types and Constants
 ******************************************************************************/

// Digital sequence residue (from easel.h)
typedef uint8_t DigitalResidue;

// Easel constants
constexpr DigitalResidue digitalResidueSentinel = 255;  // sentinel bytes 0,L+1 in a digital_sequence
constexpr DigitalResidue digitalResidueIllegal  = 254;  // input symbol is unmapped
constexpr float eslINFINITY = INFINITY;   // from C99 float.h
constexpr float eslCONST_LOG2 = 0.69314718055994529f;  // log(2.0)

/*******************************************************************************
 * 2. HMMER Constants (from p7_profile.h and related)
 ******************************************************************************/

// Transitions per node
constexpr int p7P_NTRANS = 7;

// Emissions per node (Match, Insert)
constexpr int p7P_NR = 2;

// Special states (N, E, C, J, B)
constexpr int p7P_NXSTATES = 5;
constexpr int p7P_NXTRANS = 2;  // Loop vs Move

// E-value parameters
constexpr int p7_NEVPARAM = 6;
constexpr int p7_NCUTOFFS = 6;
constexpr int p7_MAXABET = 20;

// Transition indices (0-6)
constexpr int p7P_MM = 0;  // Match->Match
constexpr int p7P_MI = 1;  // Match->Insert
constexpr int p7P_MD = 2;  // Match->Delete
constexpr int p7P_IM = 3;  // Insert->Match
constexpr int p7P_II = 4;  // Insert->Insert
constexpr int p7P_DM = 5;  // Delete->Match
constexpr int p7P_DD = 6;  // Delete->Delete

// Emission indices
constexpr int p7P_MSC = 0;  // Match Score
constexpr int p7P_ISC = 1;  // Insert Score

// Special state indices for generic DP matrix
enum p7g_xcells_e {
    p7G_E = 0,  // End
    p7G_N = 1,  // N-terminal
    p7G_J = 2,  // Join
    p7G_B = 3,  // Begin
    p7G_C = 4   // C-terminal
};
constexpr int p7G_NXCELLS = 5;

// Standard state indices for generic DP matrix
enum p7g_scells_e {
    p7G_M = 0,  // Match
    p7G_I = 1,  // Insert
    p7G_D = 2   // Delete
};
constexpr int p7G_NSCELLS = 3;

/*******************************************************************************
 * 3. HMMER Utility Macros
 ******************************************************************************/

#define ESL_MAX(a,b) (((a)>(b))?(a):(b))

#endif // MSV_FILTER_HMMER_TYPES_HPP
