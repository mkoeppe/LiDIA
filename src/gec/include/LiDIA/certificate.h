#ifndef LIDIA_CERTIFICATE_H_GUARD_
#define LIDIA_CERTIFICATE_H_GUARD_

#include        <string>

#include        "LiDIA/bigint.h"
#include        "LiDIA/LiDIA.h"
#include        "LiDIA/base_vector.h"

#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class certificate
{
private:
	base_vector <bigint>            r;

public:
	certificate();
	~certificate();
	
	void set_certificate(base_vector <bigint> vector);
	void write_certificate(std::string const& file);
	void read_certificate(std::string const& file);
	void add(certificate cert);
	void delete_all();
	base_vector <bigint> get_certification_vector();
	base_vector <bigint> get_cert_vector(int test_number);
	int get_number_of_tests();
	
};


#ifdef LIDIA_NAMESPACE
}       // end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif


#endif



