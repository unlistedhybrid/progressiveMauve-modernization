#ifndef _threadstorage_h_
#define _threadstorage_h_

// Portable thread-local storage for MUSCLE, supporting MSVC, GCC, Clang, OpenMP, single-threaded

#if defined(_MSC_VER)
    #define TLS_DECL __declspec(thread)
#elif defined(__GNUC__) || defined(__clang__)
    #define TLS_DECL __thread
#else
    #define TLS_DECL
#endif

template <typename T>
class TLS
{
public:
    TLS() {}
    explicit TLS(T t_val) { value = t_val; }

    // Use .get() everywhere for cross-platform compatibility
    T* get() const { return value; }
    void set(T* v) { value = v; }

    // Operator overloads for code using g_SeqVec->method()
    T* operator->() const { return value; }
    T& operator*()  const { return *value; }

private:
    // No copying
    TLS(const TLS&) = delete;
    TLS& operator=(const TLS&) = delete;

    TLS_DECL static T* value;
};

template <typename T>
TLS_DECL T* TLS<T>::value = nullptr;

#endif // _threadstorage_h_
