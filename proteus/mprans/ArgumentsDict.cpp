#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"

#define FORCE_IMPORT_ARRAY
#include "ArgumentsDict.h"

namespace xt
{
#if defined(__GNUC__) && !defined(__clang__)
    namespace workaround
    {
        inline void complex_allocator()
        {
           std::allocator<int> ai;
           std::allocator<double> ad;
        }
    }
#endif
}

namespace py = pybind11;

namespace proteus
{
    template <class K, class T>
    void bind_pyarray_dict(py::class_<pyarray_dict<K, T>>& cl, const std::string& name)
    {
        using map_type = pyarray_dict<K, T>;
        using key_type = typename map_type::key_type;
        using mapped_type = typename map_type::mapped_type;

        cl.def(py::init<>());

        cl.def("__setitem__",
                [](map_type& m, const key_type& k, mapped_type& v)
                {
                    auto it = m.find(k);
                    if(it != m.end())
                    {
                        it->second = std::move(v);
                    }
                    else
                    {
                        m.emplace(k, std::move(v));
                    }
                }
        );

        cl.def("__repr__",
                [name](map_type& m)
                {
                    std::ostringstream oss;
                    oss << name << '{';
                    bool f = false;
                    for(const auto& kv: m)
                    {
                        if(f)
                        {
                            oss << ", ";
                        }
                        oss << kv.first <<": " << kv.second;
                        f = true;
                    }
                    oss << '}';
                    return oss.str();
                },
               "Return the canonical representation of this map."
        );

        cl.def("__bool__",
                [](map_type& m) -> bool
                {
                    return !m.empty();
                },
                "Check wether the map is nonempty"
        );

        cl.def("__getitem__",
                [](map_type& m, const key_type& k) -> mapped_type&
                {
                    auto it = m.find(k);
                    if(it == m.end())
                    {
                        throw std::runtime_error(detail::key_error_message(k));
                    }
                    return it->second;
                }
        );
    }
}

PYBIND11_MAKE_OPAQUE(proteus::scalar_dict<double>);
PYBIND11_MAKE_OPAQUE(proteus::scalar_dict<int>);

PYBIND11_MODULE(cArgumentsDict, m)
{
    using proteus::pyarray_dict;
    using proteus::scalar_dict;
    using proteus::arguments_dict;

    xt::import_numpy();

    using dpyarray_dict = pyarray_dict<std::string, double>;
    using ipyarray_dict = pyarray_dict<std::string, int>;

    py::class_<dpyarray_dict> dad(m, "DArrayDict");
    proteus::bind_pyarray_dict(dad, "DArrayDict");

    py::class_<ipyarray_dict> iad(m, "IArrayDict");
    proteus::bind_pyarray_dict(iad, "IArrayDict");

    py::bind_map<scalar_dict<double>>(m, "DScalarDict");
    py::bind_map<scalar_dict<int>>(m, "IScalarDict");

    py::class_<arguments_dict>(m, "ArgumentsDict")
        .def(py::init<>())
        .def_readwrite("darray", &arguments_dict::m_darray)
        .def_readwrite("iarray", &arguments_dict::m_iarray)
        .def_readwrite("dscalar", &arguments_dict::m_dscalar)
        .def_readwrite("iscalar", &arguments_dict::m_iscalar);
}

