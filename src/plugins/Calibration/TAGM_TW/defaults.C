#include <sstream>
#include <iostream>

std::ofstream flog("tagm_tw_parms_defaults.out");

void defaults() {
   for (int col = 1; col <= 102; ++col) {
      flog << "0" << std::setw(15)
           << col << std::setw(15)
           << "1" << std::setw(15)
           << "-1"<< std::setw(15)
           << "0" << std::setw(15)
           << "8" << std::setw(15)
           << "0" << std::endl;
      if (col == 9 || col == 27 || col == 81 || col == 99) {
         for (int row = 1; row <= 5; ++row) {
            flog << row << std::setw(15)
                 << col << std::setw(15)
                 << "1" << std::setw(15)
                 << "-1"<< std::setw(15)
                 << "0" << std::setw(15)
                 << "8" << std::setw(15)
                 << "0" << std::endl;
         }
      }
   }
}
