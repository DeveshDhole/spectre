// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <csignal>

// Blaze 3.8.2 has a missing include that causes nvcc to fail
// (fixed in https://bitbucket.org/blaze-lib/blaze/pull-requests/63).
// We just include it here to fix the issue.
// NOLINTBEGIN
#include <blaze/util/EnableIf.h>
#include <blaze/math/traits/ReduceTrait.h>
// NOLINTEND

#ifdef __CUDA_ARCH__
// When building for Nvidia GPUs we need to disable the use of vector
// intrinsics.
#define BLAZE_USE_VECTORIZATION 0
#endif

#ifdef SPECTRE_DEBUG
#define BLAZE_THROW(EXCEPTION)           \
  struct sigaction handler {};           \
  handler.sa_handler = SIG_IGN;          \
  handler.sa_flags = 0;                  \
  sigemptyset(&handler.sa_mask);         \
  sigaction(SIGTRAP, &handler, nullptr); \
  raise(SIGTRAP);                        \
  throw EXCEPTION
#else  // SPECTRE_DEBUG
#define BLAZE_THROW(EXCEPTION)
#endif  // SPECTRE_DEBUG
