
#ifndef _UTILS_H_
#define _UTILS_H_

template <class _Fn> static void par_for(int limit, _Fn f) {
  std::vector<std::future<bool>> threads(std::thread::hardware_concurrency());
  //std::vector<std::future<bool>> threads(1);

  std::atomic<int> sequence;
  for (auto &t : threads) {
    t = std::async([&](){
      int seq;
      while ((seq = sequence++) < limit) {
        f(seq);
      }
      return true;
    });
  }

  for (auto &t : threads) {
    t.get();
  }
}

struct hex {
  char str[17];

  template <class _Ty> hex(_Ty n) {
    static_assert(sizeof(n) <= 8, "hex() needs uint64_t or less");
    std::uint64_t value = (std::uint64_t)n << (64 - sizeof(n)*8);
    for (int i = 0; i != sizeof(n)*2; ++i) {
      str[i] = "0123456789abcdef"[value >> 60];
      value <<= 4;
    }
    str[sizeof(n)*2] = 0;
  }

  operator const char *() { return str; }
};

#endif

