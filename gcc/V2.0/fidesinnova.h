#ifndef FIDESINNOVA_H
#define FIDESINNOVA_H

#include <stdint.h>
#include "polynomial.h"
#include "proofGenerator.cpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>

using namespace std;
extern "C" void store_register_instances();

class fidesinnova {
public:
  void proofGenerator();

private:
};

#endif  // FIDESINNOVA_H
