#ifndef _threadstorage_h_
#define _threadstorage_h_

#if defined(_MSC_VER)
    #define TLS_DECL __declspec(thread)
#elif defined(__GNUC__) || defined(__clang__)
    #define TLS_DECL __thread
#else
    #define TLS_DECL
#endif

// Default TLS: for pointer types
template <typename T>
class TLS
{
public:
    TLS() {}
    void set(T* v) { value = v; }
    T* get() const { return value; }
    T* operator->() const { return value; }
    T& operator*()  const { return *value; }

private:
    TLS(const TLS&) = delete;
    TLS& operator=(const TLS&) = delete;
    TLS_DECL static T* value;
};

template <typename T>
TLS_DECL T* TLS<T>::value = nullptr;

// Partial specialization for arrays
template <typename T, size_t N>
class TLS<T[N]>
{
public:
    TLS() {}
    T* get() { return value; }
    const T* get() const { return value; }
    T& operator[](size_t idx) { return value[idx]; }
    const T& operator[](size_t idx) const { return value[idx]; }
    // Assignment from array or pointer
    void set(const T* v) { for (size_t i = 0; i < N; ++i) value[i] = v[i]; }
private:
    TLS(const TLS&) = delete;
    TLS& operator=(const TLS&) = delete;
    TLS_DECL static T value[N];
};

template <typename T, size_t N>
TLS_DECL T TLS<T[N]>::value[N] = {};

#endif // _threadstorage_h_
