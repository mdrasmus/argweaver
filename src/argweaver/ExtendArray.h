/*=============================================================================

  Matt Rasmussen
  Copyright 2007-2011

  Simple vector/array class

=============================================================================*/


#ifndef SPIDIR_EXTEND_ARRAY_H
#define SPIDIR_EXTEND_ARRAY_H

#include <assert.h>
#include <algorithm>


namespace spidir {


/*=============================================================================
    equivalent of realloc for C++
*/
template <class T>
T *resize(T *array, size_t oldsize, size_t newsize)
{
    T *tmp = new T[newsize];

    // failed to get memory
    if (!tmp)
        return NULL;

    // equivalent to new
    if (oldsize == 0)
        return tmp;

    std::copy(array, array + oldsize, tmp);

    delete [] array;
    return tmp;
}


/*=============================================================================
    easy, detachable wrapper for arrays

    len      -- the length of the populated data
    capacity -- the size of the allocated array
*/
template <class ValueType>
class ExtendArray
{
public:
    typedef ValueType* ValuePtrType;
    typedef ValuePtrType* ValuePtrPtrType;

    ExtendArray(int _len=0, int _capacity=0, ValueType *_data=NULL,
                int _minsize=40) :
        data(_data),
        len(_len),
        capacity(_capacity),
        minsize(_minsize)
    {
        // make sure capacity is atleast length of data
        if (capacity < len)
            capacity = len;

        // make sure min allocate size is atleast capacity
        // useful for multiple detachment
        if (minsize < capacity)
            minsize = capacity;

        // if no pointer is given to manage, allocate our own
        if (data == NULL && capacity != 0) {
            data = new ValueType [capacity];
        }
    }

    ExtendArray(const ExtendArray &other) :
        len(other.len),
        capacity(other.capacity),
        minsize(other.minsize)
    {
        // allocate new memory
        data = new ValueType [capacity];

        // copy over data
        for (int i=0; i<len; i++)
            data[i] = other.data[i];
    }

    ~ExtendArray()
    {
        if (data)
            delete [] data;
    }

    ExtendArray &operator=(const ExtendArray &other)
    {
        ensureSize(other.len);
        len = other.len;

        // copy over data
        for (int i=0; i<len; i++)
            data[i] = other.data[i];

        return other;
    }


    bool operator==(const ExtendArray &other) const
    {
        if (len != other.len)
            return false;

        // copy over data
        for (int i=0; i<len; i++)
            if (data[i] != other.data[i])
                return false;

        return true;
    }


    ValueType *detach()
    {
        ValueType *ret = data;
        data = NULL;
        len = 0;
        capacity = 0;

        return ret;
    }

    void attach(ValueType* _data, int _len, int _capacity)
    {
        data = _data;
        len = _len;
        capacity = _capacity;
    }

    //=========================================================================
    // capacity management

    bool setCapacity(int newsize)
    {
        int oldsize = capacity;
        ValueType *ret = resize(data, oldsize, newsize);

        // failed to alloc memory
        if (!ret)
            return false;

        data = ret;
        capacity = newsize;
        return true;
    }

    bool increaseCapacity()
    {
        int newsize = capacity;
        if (newsize < minsize)
            newsize = minsize;
        newsize *= 2;
        return setCapacity(newsize);
    }

    bool ensureSize(int needed)
    {
        if (needed <= capacity)
            return true;

        int newsize = needed;
        if (newsize < minsize)
            newsize = minsize;
        while (newsize < needed)
            newsize *= 2;

        return setCapacity(newsize);
    }

    inline int get_capacity() const
    {
        return capacity;
    }


    //=========================================================================
    // data access

    inline void append(const ValueType &val)
    {
        assert(ensureSize(len + 1));
        data[len++] = val;
    }

    inline void extend(ValueType *vals, int nvals)
    {
        assert(ensureSize(len + nvals));

        for (int i=0; i<nvals; i++)
            data[len++] = vals[i];
    }

    inline ValueType pop()
    {
        assert(len > 0);
        return data[--len];
    }

    inline ValueType last()
    {
        return data[len-1];
    }

    inline void clear()
    {
        len = 0;
    }

    inline ValueType &operator[](const int i) const
    {
        return data[i];
    }

    inline ValueType *get() const
    {
        return data;
    }

    // easy access to underling data
    inline operator ValuePtrType() const
    {
        return data;
    }

    inline int size() const
    {
        return len;
    }

    inline void setSize(int _size)
    {
        len = _size;
    }

protected:
    ValueType *data;
    int len;
    int capacity;
    int minsize;
};



template <class ValueType>
class StackPointer
{
public:
    typedef ValueType* ValuePtrType;
    typedef ValueType** ValuePtrPtrType;

    StackPointer(ValueType *ptr=NULL) :
        ptr(ptr)
    {
    }

    ~StackPointer()
    {
        if (ptr)
            delete ptr;
    }

    ValueType *detach()
    {
        ValueType *ret = ptr;
        ptr = NULL;
        return ret;
    }


    ValuePtrType &get() const
    { return ptr; }

    operator ValuePtrType() const
    { return ptr; }

    ValuePtrPtrType operator &() const
    { return &ptr; }


protected:
    ValueType *ptr;
};


template <class ValueType>
class StackArray
{
public:
    typedef ValueType* ValuePtrType;
    typedef ValueType** ValuePtrPtrType;

    StackArray(ValueType *ptr=NULL) :
        ptr(ptr)
    {
    }

    ~StackArray()
    {
        if (ptr)
            delete [] ptr;
    }

    ValueType *detach()
    {
        ValueType *ret = ptr;
        ptr = NULL;
        return ret;
    }

    ValuePtrType &get() const
    { return ptr; }

    operator ValuePtrType() const
    { return ptr; }

    ValuePtrPtrType operator &()
    { return &ptr; }


protected:
    ValueType *ptr;
};



} // namespace spidir

#endif // SPIDIR_EXTEND_ARRAY_H
