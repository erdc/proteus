#ifndef ARGUMENTSDICT_H
#define ARGUMENTSDICT_H

#include <map>
#include <stdexcept>
#include <string>
#include "xtensor-python/pyarray.hpp"
#include "xtensor/xio.hpp"

namespace proteus
{


    // Special map that behaves like a standard map except for
    // - operator[] which throws instead of inserting a default
    //   value when a key is not found
    // - insert and insert_or_assign which always move the key
    //   and the value. This allows to keep pyarrays synchronized
    //   with the numpy arrays they refer to, while avoiding code
    //   duplication for array of scalars (because insert_or_assign
    //   is not available in C++14)
    template <class K, class T>
    class throwing_map : private std::map<K, T>
    {
    public:

        using base_type = std::map<K, T>;
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
        std::pair<iterator, bool> insert_or_assign(key_type&& k, mapped_type&& v);
    };
    
    template <class K, class T>
    using pyarray_dict = throwing_map<K, xt::pyarray<T>>;

    template <class K, class T>
    using scalar_dict = throwing_map<K, T>;

    struct arguments_dict
    {
        pyarray_dict<std::string, double> m_darray;
        pyarray_dict<std::string, int> m_iarray;
        scalar_dict<std::string, double> m_dscalar;
        scalar_dict<std::string, int> m_iscalar;
    };

    /*******************************
     * throwing_map implementation *
     *******************************/

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

    template <class K, class T>
    inline auto throwing_map<K, T>::operator[](const key_type& k) -> mapped_type&
    {
        auto it = this->find(k);
        if(it == end())
        {
            throw std::runtime_error(detail::key_error_message(k));
        }
        return it->second;
    }

    template <class K, class T>
    inline auto throwing_map<K, T>::insert(value_type&& v) -> std::pair<iterator, bool>
    {
        return base_type::insert(std::move(v));
    }

    template <class K, class T>
    inline auto throwing_map<K, T>::insert_or_assign(key_type&& k, mapped_type&& v) -> std::pair<iterator, bool>
    {
        auto it = base_type::find(k);
        if(it == base_type::end())
        {
            return base_type::emplace(k, std::move(v));
        }
        else
        {
            it->second = std::move(v);
            return std::make_pair(it, false);
        }
    }
}

#endif

