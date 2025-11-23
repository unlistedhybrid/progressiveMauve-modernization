#ifndef _threadstorage_h_
#define _threadstorage_h_

#define TLS_THREAD_LOCAL thread_local

// =========================
// Generic scalar or enum type
// =========================
template <typename T>
class TLS
{
public:
    TLS() {}
    explicit TLS(const T& t_val) { value = t_val; }

    TLS& operator=(const T& v) { value = v; return *this; }

    void set(const T& v) { value = v; }
    T get() const { return value; }
    T& ref() { return value; }
    const T& ref() const { return value; }

    // For pointer types (works for T = pointer)
    T* operator->() { return &value; }
    const T* operator->() const { return &value; }
    T& operator*() { return value; }
    const T& operator*() const { return value; }

private:
    TLS(const TLS&) = delete;
    TLS& operator=(const TLS&) = delete;
    static TLS_THREAD_LOCAL T value;
};

template <typename T>
TLS_THREAD_LOCAL T TLS<T>::value = T();

// =========================
// Specialization for pointers (so set/get work as expected)
// =========================
template <typename T>
class TLS<T*>
{
public:
    TLS() { value = nullptr; }
    explicit TLS(T* t_val) { value = t_val; }

    TLS& operator=(T* v) { value = v; return *this; }

    void set(T* v) { value = v; }
    T* get() const { return value; }

    T* operator->() const { return value; }
    T& operator*()  const { return *value; }

private:
    TLS(const TLS&) = delete;
    TLS& operator=(const TLS&) = delete;
    static TLS_THREAD_LOCAL T* value;
};

template <typename T>
TLS_THREAD_LOCAL T* TLS<T*>::value = nullptr;

// =========================
// Specialization for C arrays (e.g. unsigned int[256])
// =========================
template <typename T, size_t N>
class TLS<T[N]>
{
public:
    TLS() { for (size_t i = 0; i < N; ++i) value[i] = T{}; }
    explicit TLS(const T& t_val) { for (size_t i = 0; i < N; ++i) value[i] = t_val; }

    T* get() { return value; }
    const T* get() const { return value; }
    T& operator[](size_t idx) { return value[idx]; }
    const T& operator[](size_t idx) const { return value[idx]; }
    void set(const T* v) { for (size_t i = 0; i < N; ++i) value[i] = v[i]; }

private:
    TLS(const TLS&) = delete;
    TLS& operator=(const TLS&) = delete;
    static TLS_THREAD_LOCAL T value[N];
};

template <typename T, size_t N>
TLS_THREAD_LOCAL T TLS<T[N]>::value[N] = {};

#endif // _threadstorage_h_
