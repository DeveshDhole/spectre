// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <memory>
#include <optional>
#include <string>

#include "DataStructures/DataVector.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/Structure/ObjectLabel.hpp"
#include "Options/Auto.hpp"
#include "Options/Context.hpp"
#include "Options/String.hpp"
#include "Utilities/TMPL.hpp"

namespace domain::creators::time_dependent_options {
/// @{
/*!
 * \brief Read in FunctionOfTime coefficients from an H5 file and volume
 * subfile.
 *
 * \details The H5 file will only be accessed in the `retrieve_function_of_time`
 * function.
 */
struct FromVolumeFile {
  struct H5Filename {
    using type = std::string;
    static constexpr Options::String help{
        "Name of H5 file to read functions of time from."};
  };

  struct SubfileName {
    using type = std::string;
    static constexpr Options::String help{
        "Subfile that holds the volume data. Must be an h5::VolumeData "
        "subfile."};
  };

  using options = tmpl::list<H5Filename, SubfileName>;
  static constexpr Options::String help =
      "Read function of time coefficients from a volume subfile of an H5 file.";

  FromVolumeFile() = default;
  FromVolumeFile(std::string h5_filename, std::string subfile_name);

  /*!
   * \brief Searches the last observation in the volume subfile and returns
   * clones of all functions of time in \p function_of_time_names.
   *
   * \details If a value for \p time is specified, will ensure that \p time is
   * within the `domain::FunctionsOfTime::FunctionOfTime::time_bounds()` of each
   * function of time.
   */
  FunctionsOfTimeMap retrieve_function_of_time(
      const std::unordered_set<std::string>& function_of_time_names,
      const std::optional<double>& time) const;

 private:
  std::string h5_filename_;
  std::string subfile_name_;
};
/// @}
}  // namespace domain::creators::time_dependent_options
