#include "filereading.h"

#include <algorithm>
#include <iostream>

//================================
// lineProcessing
//================================

void ReadFile::lineProcessing() {
  if (line_.empty()) {
      ++currentLine_;
    } else if (line_[0] == '#') {
      ++currentLine_;
      line_ = "";
    } else {
      size_t strIter = line_.find('=');
      if (strIter == line_.size())
        throw ReadFileExcept(currentLine_, line_);
      line_ = line_.substr(strIter + 1);
      ++currentLine_;
      ++currentSignificantLine_;
    }
}

//================================
// parseFile
//================================

std::vector<std::string> ReadFile::parseFile(std::ifstream &instream) {
  std::vector<std::string> lines;
  try {
      while (true) {
          instream >> std::ws;
          if (!std::getline(instream, line_))
            break;
          lineProcessing();
          if (!line_.empty())
            lines.push_back(line_);
        }
    } catch (ReadFileExcept &e) {
      std::cerr << e.what() << std::endl;
      throw;
    }
  for (auto &x : lines) {
      x.erase(std::remove_if(x.begin(), x.end(), isspace), x.end());
    }
  return lines;
}

//================================
// ReadFileExcept ctor
//================================

ReadFile::ReadFileExcept::ReadFileExcept(int linnum, std::string message)
  : message_(std::to_string(linnum) + ": " + message) {}
