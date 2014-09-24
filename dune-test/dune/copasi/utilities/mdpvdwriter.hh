// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_COPASI_VTK_MD_PVDWRITER_HH
#define DUNE_COPASI_VTK_MD_PVDWRITER_HH

#include <vector>
#include <fstream>
#include <dune/common/deprecated.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/pdelab/instationary/onestep.hh>
#include <dune/pdelab/multidomain/vtk.hh>

#include <stdlib.h>

namespace Dune
{

    template< class MDGrid, class MultiGFS, int k >
    class MdPVDWriter
    {
        const int mNumSubDomains = k;
        MDGrid& mdgrid;
        MultiGFS& multigfs;
        std::string basename;
        PDELab::FilenameHelper fn;
        std::string path;
        std::vector<double> timesteps;
        Dune::VTK::OutputType outputtype;
        Dune::VTK::DataMode datamode;
        unsigned int offset;
        unsigned int currentOffset;

    public:

        MdPVDWriter(MDGrid & mdgrid_, MultiGFS& multigfs_, std::string basename_,
                  Dune::VTK::DataMode datamode_ = Dune::VTK::conforming,
                  Dune::VTK::OutputType outputtype_ = Dune::VTK::appendedraw,
                  std::string path_="vtk", unsigned int offset_=0)
            : mdgrid(mdgrid_)
            , multigfs(multigfs_)
            , basename(basename_)
            , fn(basename_,offset_)
            , path(path_)
            , timesteps()
            , outputtype(outputtype_)
            , datamode(datamode_)
            , offset(offset_)
            , currentOffset(offset_)
        {}

        template<class V> void write(double time, const V& unew)
        {
            /* remember current time step */
            timesteps.push_back(time);

            for (int index = 0; index < mNumSubDomains; ++index)
            {

                std::stringstream str;
                str << "_domain" << index;

                PDELab::FilenameHelper currentFn(basename + str.str(), currentOffset);

                typedef typename MDGrid::SubDomainGrid SDGrid;
                const SDGrid& sdgrid = mdgrid.subDomain(index);
                typename SDGrid::LeafGridView sdgv = sdgrid.leafGridView();

                Dune::VTKWriter<typename SDGrid::LeafGridView> vtkwriter(sdgv,datamode);
                Dune::PDELab::MultiDomain::addSolutionToVTKWriter(
                  vtkwriter, multigfs, unew,
                  Dune::PDELab::MultiDomain::subdomain_predicate<typename MDGrid::SubDomainIndex>
                            (sdgrid.domain())
                );
                /* write VTK file */
                vtkwriter.pwrite(currentFn.getName(),path,"",outputtype);

                /* write pvd file */
                std::string pvdname = basename + str.str() + ".pvd";
                std::ofstream pvd(pvdname.c_str());
                assert(pvd.is_open());
                pvd << std::fixed;
                pvd << "<?xml version=\"1.0\"?>\n"
                    << "<VTKFile type=\"Collection\" version=\"0.1\">\n"
                    << "<Collection>\n";
                PDELab::FilenameHelper fnloop(basename + str.str(),offset);
                for (unsigned int i=0; i<timesteps.size(); i++)
                {
                    std::string fname = vtkwriter.getParallelHeaderName(fnloop.getName(), path, sdgv.comm().size());
                    pvd << "  <DataSet timestep=\"" << timesteps[i]
                           << "\" file=\"" << fname << "\"/>\n";
                    fnloop.increment();
                }
                pvd << "</Collection>\n"
                    << "</VTKFile>\n";
                pvd.close();

            }

            /* increment counter */
            fn.increment();
            ++currentOffset;
        }
    };

} // end namespace Dune

#endif // DUNE_COPASI_VTK_MD_PVDWRITER_HH
