#ifndef REWEIGHTINGS_HANDLER_H
#define REWEIGHTINGS_HANDLER_H

#include "matrix.h"
#include <vector>

 // *********************************************************************************
 // *                                                                               *
 // *  "ReweightingsHandler" provides access to reweighting factors stored in files.*
 // *  Currently the files supported are those output from openQCD. More info about *
 // *  openQCD can be found here: http://luscher.web.cern.ch/luscher/openQCD/.      *
 // *  The particular format of the file that contains the reweighting factors must *
 // *  be specified with a tag <Format> (see required XML below).                   *
 // *                                                                               *
 // *  Reweighting factors are output in binary by the program 'ms1', which is a    *
 // *  part of openQCD. The structure of these binary files, along with information *
 // *  on using 'ms1' can be found in 'main/README.ms1' from the openQCD source.    *
 // *  If you are using these files, you must pass 'OPENQCD' to the <Format> tag.   *
 // *  Note that versions 1.2 or older of openQCD did not specify the individual    *
 // *  Hasenbusch factors, and if the reweighting factors were computed with one of *
 // *  these older versions of openQCD, you must pass 'OPENQCD_12' to the <Format>  *
 // *  tag.                                                                         *
 // *                                                                               *
 // *  openQCD also provides a program for converting these binary files to ASCII,  *
 // *  which can be found here 'devel/nompi/main/read2.c' from the openQCD source.  *
 // *  You can also use these ASCII files, in which case you must pass 'ASCII' to   *
 // *  the <Format> tag. The files containing the reweighting factors can be        *
 // *  specified with the <MCObservables> tag that occurs within the <Initialize>   *
 // *  tag of a sigmond input file. The XML has the required form:                  *
 // *                                                                               *
 // *       <ReweightingData>                                                       *
 // *          <Format>OPENQCD</Format> (or OPENQCD_12 or ASCII)                    *
 // *          <FileName>...</FileName> (order by config number)                    *
 // *             ....                                                              *
 // *       </ReweightingData>                                                      *
 // *                                                                               *
 // *  Also note that the ordering of the <FileName> tags is important. Generally,  *
 // *  you will only have multiple files if there is more than one chain for the    *
 // *  ensemble being used. For example, if you have two chains r1 and r2 with N1   *
 // *  and N2 configs, respectively, then specifying the reweighting factors file   *
 // *  associated with r1 first will associate the N1 reweighting factors of r1     *
 // *  with the first N1 configs. Then, specificying the reweighting factors file   *
 // *  associated with r2 next will associate the N2 reweighting factors of r2 with *
 // *  the next N2 configs.                                                         *
 // *                                                                               *
 // *  There is one other things that should be noted. In principle, there should   *
 // *  be no difference between using the binary file or the ascii file that        *
 // *  openQCD produces. However, there is one difference that is implemented in    *
 // *  the 'devel/nompi/main/read2.c' file, namely the use of 'lnm' on line 340 of  *
 // *  that file. In the interest of full disclosure, I don't currently fully       *
 // *  understand the motivations for this choice made by Luescher. But, according  *
 // *  to Daniel Mohler, the use of 'lnm' is not necessary and he does not use it   *
 // *  himself. I would like to investigate this further at some point.             *
 // *                                                                               *
 // *  This class provides access to the reweighting factors for the class          *
 // *  MCObsGetHandler through calls to 'getData'.                                  *
 // *                                                                               *
 // *  A general note about reweighting in sigmond:                                 *
 // *  Reweighting is applied to simple binned observables, and reweighting         *
 // *  requires access to the unrebinned data. Thus, you can only reweight          *
 // *  quantities that are stored in files on disk in an unrebinned format. You     *
 // *  cannot reweight bins stored in memory, because they will in general be       *
 // *  rebinned. Of course, a check could have been made to see if any rebinning    *
 // *  has been done, but it didn't seem worth it to implement this limited use     *
 // *  case.                                                                        *
 // *                                                                               *
 // *********************************************************************************

class ReweightingsHandler
{
    std::string file_format;
    std::vector<std::string> file_names;
    uint number_of_meas;

    int32_t nrw;
    std::vector<int32_t> nfct;
    std::vector<int32_t> nsrc;

    std::vector<std::vector<std::vector<std::vector<double> > > > sqn;
    std::vector<std::vector<std::vector<std::vector<double> > > > lnr;

    std::vector<std::vector<double> > rw;

    RVector rw_prod;

 public:

    ReweightingsHandler(const std::string in_file_format,
                        const std::vector<std::string>& in_file_names,
                        const uint in_number_of_meas)
      : file_format(in_file_format), file_names(in_file_names), 
        number_of_meas(in_number_of_meas), nrw(0) {}

    void getData(RVector& result);

    bool queryData();

 private:

    void read_data();

    void read_openqcd(const std::string file_name);

    void read_openqcd12(const std::string file_name);

    void read_ascii(const std::string file_name);

    void read_openqcd_binary(const std::string file_name, bool v12=false);
};

#endif
