typedef struct hitTree_s {
  int mark;
  struct hitTree_s* left;
  struct hitTree_s* right;
  void* this;
} hitTree_t;

void** getTwig(hitTree_t** tree, int mark);

void gmtod_(float xin[3], float xout[3], const int* xord);

