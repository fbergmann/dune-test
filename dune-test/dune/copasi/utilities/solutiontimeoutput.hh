#ifndef DUNE_TEST_SOLUTIONTIMEOUTPUT_H
#define DUNE_TEST_SOLUTIONTIMEOUTPUT_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <fstream>
#include <dune/geometry/referenceelements.hh>
#include <dune/pdelab/common/geometrywrapper.hh>


/** \brief Class to write a solution in time in 1 cell to file.
 *
 * It works for arbitrary number of DiscreteGridFunctions
 *
 * \tparam RF  C++ type of the floating point parameters
 * \tparam GFS GridFunctionSpace
 */
template<typename RF, typename GV>
class GnuplotSolutionTimeOutput {

  // get some types
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
  typedef typename GV::Traits::template Codim<0>::Entity Element;
  typedef typename Dune::PDELab::ElementGeometry<Element> EG;
  typedef std::size_t size_type;
  typedef typename GV::Grid::LeafIndexSet IndexSet;
  typedef typename IndexSet::IndexType IndexType;

  const GV &gv;
  std::vector<ElementIterator> egstorage;
  std::ofstream s;
  bool enabled;

  enum {dim = GV::dimension};

public:
  //! Construct
  /**
   * \todo            It does NOT work in parallel.
   * \param fname     Name of the file to write to.
   * \param gv_       Grid View used for the vectors passed in.

   * \note The filename passed should be the same on all ranks.  Only rank
   *       0 will actually write to the file.
   */
  GnuplotSolutionTimeOutput(const std::string &fname, const GV &gv_, const std::string &comment = " ") :
    gv(gv_)
  {

    enabled = fname.size() > 0;
    if(enabled) {

        if(gv.comm().rank() == 0) {
            s.exceptions(std::ios_base::badbit | std::ios_base::eofbit
                         | std::ios_base::failbit);

            const char * cname = fname.c_str();

            s.open(cname);
            s << std::setprecision(14) << std::scientific;
            if(comment.size() > 0)
              s << "# " << comment << "\n";
          }
      }

    std::cout << "size of vector is " << egstorage.size() << std::endl;
    // loop once over the grid and find elements with given x_coord
    for (ElementIterator it = gv.template begin<0>();
         it!=gv.template end<0>(); ++it)
      {
        EG eg(*it);

        // we will only add first element to the output storage
        // it can be modified and add to storage arbitrary element of your choice,
        // but only ONE!
        if (egstorage.size()<1)
          egstorage.push_back(it);

        //     // cell geometry
        //     const Dune::FieldVector<RF,dim>&
        //       cell_center_local = Dune::GenericReferenceElements<RF,dim>::general(eg.geometry().type()).position(0,0);
        //
        //     Dune::FieldVector<RF, dim>
        //       cell_center_global = eg.geometry().global(cell_center_local);
        //
        //     size_type nr_corners = eg.geometry().corners();
        //     size_type nr_corners_half = nr_corners/2;
        //
        //     size_type nrx = 0;
        //     size_type nry = 0;
        //
        //     for (unsigned int i=0; i<nr_corners; i++)
        //       {
        //         Dune::FieldVector<RF,dim> corner = eg.geometry().corner(i);
        //         if (corner[0] < x_coord)
        //           ++nrx;
        //         if (corner[0] > x_coord)
        //           ++nry;
        //       }
        //
        //     // stores appropriate element iterators
        //     if (nrx == nry  && nry== nr_corners_half)
        //       egstorage.push_back(it);
      }// end it

    if (egstorage.size()!=1)
      DUNE_THROW(Dune::Exception,"GnuplotSolutionTimeOutput: in egstorage should be exactly 1 Element!");

  }

  //! Write a record in the output file (Element index)
  /**
   * \param t        time
   * \param ... dgf  discrete grid function ellipses

   * example: output(100, c1dgf, c2dgf, c3dgf)
   */
  template<class Time, typename... DGF>
  void write(Time t, const DGF&... dgf) {

    // loop over all elements from egstorage
    for(typename std::vector<ElementIterator>::iterator pit=egstorage.begin();pit!=egstorage.end();++pit)
      {
        EG eg(*(*pit));

        s << t << "\t";
        DGFoutput(eg,dgf...);
        s << "\n";
      }
  }

private:

  //! output from DGF (only rank 0), ellipses
  template<typename EG, typename HDGF, typename... DGF>
  void DGFoutput(const EG& eg, const HDGF& hdgf, const DGF&... dgf) {

    typename HDGF::Traits::RangeType value;
    // cell geometry
    const Dune::FieldVector<RF,dim>&
        cell_center_local = Dune::GenericReferenceElements<RF,dim>::general(eg.geometry().type()).position(0,0);
    hdgf.evaluate(eg.entity(),cell_center_local,value);
    s << value << "\t";

    DGFoutput(eg,dgf...);
  }

  //! output from DGF
  template<typename EG, typename HDGF>
  void DGFoutput(const EG& eg, const HDGF& hdgf) {
    typename HDGF::Traits::RangeType value;
    // cell geometry
    const Dune::FieldVector<RF,dim>&
        cell_center_local = Dune::GenericReferenceElements<RF,dim>::general(eg.geometry().type()).position(0,0);
    hdgf.evaluate(eg.entity(),cell_center_local,value);
    s << value << "\t";
  }

};

#endif //DUNE_TEST_SOLUTIONTIMEOUTPUT_H


// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
