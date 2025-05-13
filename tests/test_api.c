/**
 * @file test_api.c
 * @brief Bootstraps the test framework for the CCA library.
 * @author - SCEC
 * @version 1.0
 *
 * Tests the CCA library by loading it and executing the code as
 * UCVM would do it.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "cca.h"

/**
 * Initializes and runs the test program. Tests link against the
 * static version of the library to prevent any dynamic loading
 * issues.
 *
 * @param argc The number of arguments.
 * @param argv The argument strings.
 * @return A zero value indicating success.
 */
int main(int argc, const char* argv[]) {

	// Declare the structures.
	cca_point_t pt;
	cca_properties_t ret;

        // Initialize the model.
        char *envstr=getenv("UCVM_INSTALL_PATH");
        if(envstr != NULL) {
           assert(cca_init(envstr, "cca") == 0);
           } else {
             assert(cca_init("..", "cca") == 0);
        }

	printf("Loaded the model successfully.\n");

	// Query a point.
	pt.longitude = -120;
	pt.latitude = 37.2;
	pt.depth = 0;

	cca_query(&pt, &ret, 1);

	assert(ret.vs > 0);
	assert(ret.vp > 0);
	assert(ret.rho > 0);

	printf("Query was successful.\n");

	// Close the model.
	assert(cca_finalize() == 0);

	printf("Model closed successfully.\n");

	printf("\nALL CCA TESTS PASSED\n");

	return 0;
}
