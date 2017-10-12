#ifndef SRC_FILEREADING_H_
#define SRC_FILEREADING_H_

#include <boost/noncopyable.hpp>

#include <exception>
#include <fstream>
#include <string>
#include <vector>

//================================
// ReadFile
//================================

class ReadFile : private boost::noncopyable {
  size_t currentLine_,
         currentSignificantLine_;
  std::string line_;

public:
  std::vector<std::string> parseFile(std::ifstream &instream);

public:
  class ReadFileExcept;

private:
  void lineProcessing();
};

//================================
// ReadFileExcept
//================================

class ReadFile::ReadFileExcept : public std::exception {
  std::string message_;

public:
  ReadFileExcept(int linnum, std::string message);

  virtual const char *what() const noexcept {
    return message_.c_str();
  }

  ~ReadFileExcept() noexcept {}
};
#endif  // SRC_FILEREADING_H_
