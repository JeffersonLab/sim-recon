typedef struct hitTree_s {
  int mark;
  struct hitTree_s* left;
  struct hitTree_s* right;
  void* this;
} binTree_t;

void** getTwig(binTree_t** tree, int mark);
void* pickTwig(binTree_t** tree);
