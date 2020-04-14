#ifndef ARGUMENTSDICT_H
#define ARGUMENTSDICT_H

#include <map>
#include <stdexcept>
#include <string>
#include "xtensor-python/pyarray.hpp"
#include "xtensor/xio.hpp"

namespace proteus
{
    namespace detail
    {
        template <class T>
        inline std::string key_error_messagei(const T& key)
        {
            return key_error_message(std::to_string(key));
        }

        inline std::string key_error_message(const std::string& key)
        {
            return std::string("key ") + key + std::string(" not found");
        }
    }

    // Special map for pyarrays: the arrays are moved
    // into the map instead of being copied. This allows
    // to keep them synchronized with the numpy arrays
    // they refer to.
    template <class K, class T>
    class pyarray_dict : private std::map<K, xt::pyarray<T>>
    {
    public:

        using base_type = std::map<K, xt::pyarray<T>>;
        using key_type = typename base_type::key_type;
        using mapped_type = typename base_type::mapped_type;
        using value_type = typename base_type::value_type;
        using reference = typename base_type::reference;
        using const_reference = typename base_type::const_reference;
        using pointer = typename base_type::pointer;
        using const_pointer = typename base_type::const_pointer;
        using size_type = typename base_type::size_type;
        using difference_type = typename base_type::difference_type;
        using iterator = typename base_type::iterator;
        using const_iterator = typename base_type::const_iterator;

        using base_type::empty;
        using base_type::find;
        using base_type::emplace;

        using base_type::begin;
        using base_type::end;
        using base_type::cbegin;
        using base_type::cend;

        // Throws if k is not in the map
        mapped_type& operator[](const key_type& k);
    
        std::pair<iterator, bool> insert(value_type&&);
    };
    
    template <class K, class T>
    inline auto pyarray_dict<K, T>::operator[](const key_type& k) -> mapped_type&
    {
        auto it = this->find(k);
        if(it == end())
        {
            throw std::runtime_error(detail::key_error_message(k));
        }
        return it->second;
    }

    template <class K, class T>
    inline auto pyarray_dict<K, T>::insert(value_type&& v) -> std::pair<iterator, bool>
    {
        return base_type::insert(std::move(v));
    }

    template <class T>
    using scalar_dict = std::map<std::string, T>;

    struct arguments_dict
    {
        pyarray_dict<std::string, double> m_darray;
        pyarray_dict<std::string, int> m_iarray;
        scalar_dict<double> m_dscalar;
        scalar_dict<int> m_iscalar;
    };
}

#endif

