#include <map>
#include <iostream>

extern "C" {
#include "gid_map.h"
}

static std::map <int, int> gid2id;

extern "C" {

  int gidGetId(int gid) {
    int id = gid2id[gid];
    if (id == 0) {
      id = -1;
      gid2id[gid] = id;
    }
    return id;
  }

  void gidSet(int gid, int id) {
    gid2id[gid] = id;
    return;
  }

  void gidClear() {
    gid2id.clear();
    return;
  }

  void gidclear_() {
    gidClear();
    return;
  }

}

