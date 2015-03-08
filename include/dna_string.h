
#ifndef _DNA_STRING_H_
#define _DNA_STRING_H_

#include <string>

namespace std {
  class dna_base {
    dna_base() : value(0) {}
    ~dna_base() = default;
    explicit dna_base(char x) : value(x == 'C' ? 1 : x == 'G' ? 2 : x == 'T' ? 3 : 0) {}
    explicit dna_base(unsigned char x) : value(x) {}
    operator char() const { return "ACGT"[value]; }
    dna_base& operator=(const dna_base x) { value = x.value; return *this; }
    dna_base& operator=(const dna_base& x) { value = x.value; return *this; }
    void flip() { value ^= 3; }
  private:
    unsigned char value : 2;
  };

  template <class Allocator> class vector<dna_base, Allocator> {
  public:
    // types:
    typedef dna_base const_reference;
    typedef implementation-defined iterator; // see 23.2
    typedef implementation-defined const_iterator; // see 23.2
    typedef implementation-defined size_type; // see 23.2
    typedef implementation-defined difference_type;// see 23.2
    typedef dna_base value_type;
    typedef Allocator allocator_type;
    typedef implementation-defined pointer;
    typedef implementation-defined const_pointer;
    typedef std::reverse_iterator<iterator> reverse_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

    // bit reference:
    class reference {
      friend class vector;
      reference() = default;
      reference(size_t i, _ChunkType *b) :
        _Index(i), _Begin(b) {}
    public:
      ~reference() = default;
      operator dna_base() const {
        unsigned char v = (unsigned char)(
          _Begin[_Index / _NumBasesPerChunk] >>
          (_Index % _NumBasesPerChunk)*2
        ) & 3;
        return dna_base(v);
      }
      reference& operator=(const dna_base x) {
        _ChunkType mask =
          (_ChunkType)3 << _NumBasesPerChunk*2 >>
          (_Index % _NumBasesPerChunk)*2
        ;
        _ChunkType val =
          (_ChunkType)(unsigned char)x << _NumBasesPerChunk*2 >>
          (_Index % _NumBasesPerChunk)*2
        ;
        _Begin[_Index / _NumBasesPerChunk] &= ~mask;
        _Begin[_Index / _NumBasesPerChunk] |= val;
      }
      reference& operator=(const reference& x) {
        *this = (dna_base)x;
      }
      void flip() {
        _ChunkType mask =
          (_ChunkType)3 << _NumBasesPerChunk*2 >>
          (_Index % _NumBasesPerChunk)*2
        ;
        _Begin[_Index / _NumBasesPerChunk] ^= mask;
      }
    private:
      size_t _Index;
      _ChunkType *_Begin;
    };

    // construct/copy/destroy:
    explicit vector(const Allocator& = Allocator());
    explicit vector(size_type n, const dna_base& value = dna_base(),

    const Allocator& = Allocator());
    template <class InputIterator>
    vector(InputIterator first, InputIterator last,
    const Allocator& = Allocator());
    vector(const vector<dna_base,Allocator>& x);
    vector(vector<dna_base,Allocator>&& x);
    vector(const vector&, const Allocator&);
    vector(vector&&, const Allocator&);
    vector(initializer_list<dna_base>, const Allocator& = Allocator()));
    ~vector();
    vector<dna_base,Allocator>& operator=(const vector<dna_base,Allocator>& x);
    vector<dna_base,Allocator>& operator=(vector<dna_base,Allocator>&& x);
    vector operator=(initializer_list<dna_base>);
    template <class InputIterator>
    void assign(InputIterator first, InputIterator last);
    void assign(size_type n, const dna_base& t);
    void assign(initializer_list<dna_base>;
    allocator_type get_allocator() const;
    // iterators:
    iterator begin();
    const_iterator begin() const;
    iterator end();
    const_iterator end() const;
    reverse_iterator rbegin();
    const_reverse_iterator rbegin() const;
    reverse_iterator rend();
    const_reverse_iterator rend() const;
    const_iterator cbegin() const;
    const_iterator cend() const;
    const_reverse_iterator crbegin() const;
    const_reverse_iterator crend() const;
    // capacity:
    size_type size() const;
    size_type max_size() const;
    void resize(size_type sz, dna_base c = false);
    size_type capacity() const;
    bool empty() const;
    void reserve(size_type n);
    void shrink_to_fit();
    // element access:
    reference operator[](size_type n);
    const_reference operator[](size_type n) const;
    const_reference at(size_type n) const;
    reference at(size_type n);
    reference front();
    const_reference front() const;
    reference back();
    const_reference back() const;
    // modifiers:
    void push_back(const dna_base& x);
    void pop_back();
    iterator insert(const_iterator position, const dna_base& x);
    iterator insert (const_iterator position, size_type n, const dna_base& x);
    template <class InputIterator>
    iterator insert(const_iterator position,
    InputIterator first, InputIterator last);
    iterator insert(const_iterator position, initializer_list<dna_base> il);
    iterator erase(const_iterator position);
    iterator erase(const_iterator first, const_iterator last);
    void swap(vector<dna_base,Allocator>&);
    static void swap(reference x, reference y);
    void flip(); // flips all bits
    void clear();

  private:
    typedef unsigned long long _ChunkType;
    static const size_t _NumBasesPerChunk = sizeof(_ChunkType)/2;
    _ChunkType *_Begin;
    _ChunkType *_End;
    _ChunkType *_Capacity;
  };
}


#endif

