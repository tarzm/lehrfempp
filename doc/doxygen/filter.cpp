// Author: Huang Liaowang
// Date: January 2020
// This is part of the LehrFEM++ code suite

/*
This program replaces all the @lref{<label>} statements in a C++ files
with a label-specific string and the corresponding number from .aux file. Here
we assume the line corresponding to this label in .aux file looks like

\newlabel{<label>@cref}{{[<Label>][x][xxx]<number>}{[x][xx][x]xxx}}

We would like to replace @lref{<label>} with <Label> <number>, which means that
label is extracted from the first square brackets and the number comes from the
content after two brackets of <Label>.

This code is heavily optimized for speed. In particular it will create a cache
file which contains a lookup-table of the .aux file so that we can map <label>
-> <Label> <number>
*/

/*
The ultimate goal is to have references to LaTeX documents in documentation
generated by Doxygen, and this script is to do the preprocessing before the
Doxygen parse the source file.
*/

#include <boost/iostreams/device/mapped_file.hpp>
#include <fstream>
#include <iostream>
#include <map>
#include <optional>
#include <regex>
#include <string>
#include <string_view>
#include <unordered_map>

#include <boost/fusion/adapted.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/karma_string.hpp>
#include <boost/spirit/include/qi.hpp>

class Filter {
 private:
  unordered_map<std::string, std::string> label_map_;
  std::vector<std::pair<std::string, std::string>> aux_table_;

  std::string aux_file_;
  std::string file_;

 public:
  Filter(std::string file, std::string aux_file)
      : file_(std::move(file)), aux_file_(std::move(aux_file)) {
    label_map_ = {{"equation", "Equation"},
                  {"par", "Paragraph"},
                  {"chapter", "Chapter"},
                  {"sec", "Section"},
                  {"subsection", "Subsection"},
                  {"figure", "Figure"},
                  {"code", "Code"},
                  {"remark", "Remark"},
                  {"subsubsection", "Subsection"},
                  {"example", "Example"}};
  }

  /**
   * @brief Reads a file from disk into a std::string. If the file doesn't
   * exist or cannot be opened, an empty optional is returned.
   * @param file_name The path to the file that should be opened.
   */
  std::optional<std::string> ReadFile(const std::string& file_name);

  /**
   * @brief Load the .aux file and create the `aux_table_`
   *
   * The line in .aux file we are interested is assumed to look like:
   * `\newlabel{<label>@cref}{{[<Label>][x][xxx]<number>}{[x][xx][x]xxx}}`
   *
   *
   * If there is a cache file present, the hash of the `.aux` file is compared
   * to the stored hash in the cache file. If they agree, the lookup table is
   * loaded from cache. Otherwise the lookup table (`aux_table_`) is created
   * from the `.aux` file.
   */
  void LoadAuxTable();
  /**
   * @brief Goes through the source file passed by doxygen (second commandline
   * argument) and replace @lref{} with the corresponding value in the lookup
   * table (`aux_table_`). The result is output to the console (`std::cout`)
   */
  void ReplaceLref();
};

std::optional<std::string> Filter::ReadFile(const std::string& file_name) {
  ifstream t(file_name, ios::binary);  // open the file
  if (!t.is_open()) {
    return {};
  }
  t.seekg(0, std::ios::end);
  std::size_t aux_size = t.tellg();
  std::string result;
  result.resize(aux_size);
  t.seekg(0);
  t.read(&result[0], aux_size);

  return result;
}

void Filter::LoadAuxTable() {
  namespace qi = boost::spirit::qi;
  namespace karma = boost::spirit::karma;

  // load auxilliary file and compute it's hash
  auto aux_string = ReadFile(aux_file_);
  if (!aux_string.has_value()) {
    std::cerr << "filter.cpp : Could not open " << aux_file_ << std::endl
              << std::flush;
    std::exit(1);
  }
  auto aux_hash = std::hash<std::string>()(aux_string.value());

  // load cache file and retrieve hash
  std::string cache_filename = "filter.cache";
  bool rebuild_cache = false;  // should we save our cache at the end?
  auto cache_string = ReadFile(cache_filename);
  if (cache_string.has_value()) {
    std::size_t cache_hash;
    std::size_t num_entries;
    auto begin = cache_string.value().cbegin();
    auto end = cache_string.value().cend();
    if (!qi::phrase_parse(begin, end, qi::ulong_long >> qi::ulong_long, qi::eol,
                          cache_hash, num_entries)) {
      std::cerr << "filter.cpp : Cannot get hash from " << cache_filename
                << " -> rebuilding cache" << std::endl;
      rebuild_cache = true;
    } else {
      if (cache_hash != aux_hash) {
        // if hash stored in cache file doesn't agree with our computed hash
        // of the aux file, rebuild cache
        rebuild_cache = true;
        std::cerr << "filter.cpp : Hashes don't agree -> Rebuild Cache"
                  << std::endl;
      } else {
        // load the rest of the cache
        aux_table_.reserve(num_entries);
        if (!qi::phrase_parse(begin, end,
                              *(qi::lexeme[+(qi::char_ - qi::eol)] >>
                                qi::lexeme[+(qi::char_ - qi::eol)]),
                              qi::eol, aux_table_)) {
          std::cerr
              << "filter.cpp : Cannot load cache entries -> rebuilding cache";
          rebuild_cache = true;
        }
      }
    }
  } else {
    rebuild_cache = true;
  }

  if (rebuild_cache) {
    // build cache:
    auto begin = aux_string.value().cbegin();
    smatch matches;
    regex label_pattern(
        R"(\\newlabel\{(.*)@cref\}\{\{\[([^,\]]*)\](\[[^\}]*\])*([0-9\.]*)\})",
        std::regex::optimize);
    while (regex_search(begin, aux_string.value().cend(), matches,
                        label_pattern)) {
      if (matches.ready() && !matches.empty()) {
        std::string label = matches.str(1);
        // NOLINTNEXTLINE
        std::string Label = label_map_.count(matches.str(2)) > 0U
                                ? label_map_[matches.str(2)]
                                : matches.str(2);
        std::string number = matches.str(4);

        if (!matches.str(1).empty()) {
          // convert the <Label> into human readable form:
          // aux_table_[matches.str(1)] = Label + " " + number;
          // NOLINTNEXTLINE
          aux_table_.emplace_back(matches.str(1), Label + " " + number);
        }
        begin += matches.prefix().length() + matches[0].length();
      }
    }
    std::sort(aux_table_.begin(), aux_table_.end(),
              [](auto& a, auto& b) { return a.first < b.first; });

    // save cache
    std::ofstream cache_out(cache_filename, std::ios_base::out |
                                                std::ios_base::binary |
                                                std::ios_base::trunc);
    karma::ostream_iterator<char> outit(cache_out);

    if (!karma::generate_delimited(
            outit,
            karma::ulong_long << karma::ulong_long
                              << *(karma::string << karma::string),
            karma::eol, aux_hash, aux_table_.size(), aux_table_)) {
      std::cerr << "filter.cpp : Error while writing to cache!" << std::endl;
      std::exit(1);
    }
  }
}

void Filter::ReplaceLref() {
  // disable sync with C streams -> faster output to command line.
  std::ios::sync_with_stdio(false);

  // Read file into memory
  auto file_string = ReadFile(file_);
  if (!file_string.has_value()) {
    std::cerr << "filter.cpp : Could not open " << file_ << std::endl
              << std::flush;
    std::exit(1);
  }

  // construct result with regex.
  std::string result;
  result.reserve(file_string.value().size() + 100);
  regex lref("@lref\\{(.*?)\\}", regex_constants::optimize);
  smatch match;
  auto begin = file_string.value().cbegin();
  while (regex_search(begin, file_string.value().cend(), match, lref)) {
    result += match.prefix().str();
    if (aux_table_.empty()) {
      // build aux lookup table only if we need it.
      LoadAuxTable();
    }
    auto search_result =
        std::lower_bound(aux_table_.begin(), aux_table_.end(),
                         std::pair<std::string, std::string>(match[1], ""),
                         [](const std::pair<std::string, std::string>& a,
                            const std::pair<std::string, std::string>& b) {
                           return a.first < b.first;
                         });
    if (search_result == aux_table_.end()) {
      cerr << "filter.cpp : Cannot find a label for " << match[1] << endl;
      result += match[0];
    } else {
      result += search_result->second;
    }
    begin += match.prefix().length() + match[0].length();
  }
  result.append(begin, file_string.value().cend());

  std::cout << result;
}

int main(int argc, char* argv[]) {
  // takes a file as input and print the content after parsing.
  // the output is cout which can be understood natively by doxygen.
  // Process arguments
  std::string aux_file;
  std::string file;
  switch (argc) {
    case 1: {
      std::cout << "Usage: " << argv[0] << " [aux file] <source file>" << endl;
      return -1;
    }
    case 2: {
      aux_file = "NPDEFLrefs.aux";
      file = argv[1];
      break;
    }
    default: {
      aux_file = argv[1];
      file = argv[2];
      break;
    }
  }

  // Perform filtering
  Filter f(file, aux_file);
  f.ReplaceLref();

  return 0;
}
