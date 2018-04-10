/**
 * @file cca.c
 * @brief Main file for CCA library.
 * @author David Gill - SCEC <davidgil@usc.edu>
 * @version 1.0
 *
 * @section DESCRIPTION
 *
 * Delivers the prototype CCA model which consists of En-Jui Lee's full 3D
 * tomographic results for central California.
 *
 */

#include "cca.h"

/**
 * Initializes the CCA plugin model within the UCVM framework. In order to initialize
 * the model, we must provide the UCVM install path and optionally a place in memory
 * where the model already exists.
 *
 * @param dir The directory in which UCVM has been installed.
 * @param label A unique identifier for the velocity model.
 * @return Success or failure, if initialization was successful.
 */
int cca_init(const char *dir, const char *label) {
	int tempVal = 0;
	char configbuf[512];
	double north_height_m = 0, east_width_m = 0, rotation_angle = 0;

	// Initialize variables.
	cca_configuration = calloc(1, sizeof(cca_configuration_t));
	cca_velocity_model = calloc(1, sizeof(cca_model_t));
        cca_vs30_map = calloc(1, sizeof(cca_vs30_map_config_t));

	// Configuration file location.
	sprintf(configbuf, "%s/model/%s/data/config", dir, label);

        // Set up model directories.
        sprintf(vs30_etree_file, "%s/model/ucvm/ucvm.e", dir);

	// Read the cca_configuration file.
	if (read_cca_configuration(configbuf, cca_configuration) != SUCCESS)
		return FAIL;

	// Set up the iteration directory.
	sprintf(cca_iteration_directory, "%s/model/%s/data/%s/", dir, label, cca_configuration->model_dir);

	// Can we allocate the model, or parts of it, to memory. If so, we do.
	tempVal = try_reading_model(cca_velocity_model);

	if (tempVal == SUCCESS) {
		fprintf(stderr, "WARNING: Could not load model into memory. Reading the model from the\n");
		fprintf(stderr, "hard disk may result in slow performance.");
	} else if (tempVal == FAIL) {
		print_error("No model file was found to read from.");
		return FAIL;
	}

        if (read_vs30_map(vs30_etree_file, cca_vs30_map) != SUCCESS) {
                print_error("Could not read the Vs30 map data from UCVM.");
                return FAIL;
        }

	// We need to convert the point from lat, lon to UTM, let's set it up.
	if (!(cca_latlon = pj_init_plus("+proj=latlong +datum=WGS84"))) {
		print_error("Could not set up latitude and longitude projection.");
		return FAIL;
	}
	if (!(cca_utm = pj_init_plus("+proj=utm +zone=11 +ellps=clrk66 +datum=NAD27 +units=m +no_defs"))) {
		print_error("Could not set up UTM projection.");
		return FAIL;
	}
        if (!(cca_aeqd = pj_init_plus(cca_vs30_map->projection))) {
                print_error("Could not set up AEQD projection.");
                return FAIL;
        }

	// In order to simplify our calculations in the query, we want to rotate the box so that the bottom-left
	// corner is at (0m,0m). Our box's height is total_height_m and total_width_m. We then rotate the
	// point so that is is somewhere between (0,0) and (total_width_m, total_height_m). How far along
	// the X and Y axis determines which grid points we use for the interpolation routine.

	// Calculate the rotation angle of the box.
	north_height_m = cca_configuration->top_left_corner_n - cca_configuration->bottom_left_corner_n;
	east_width_m = cca_configuration->top_left_corner_e - cca_configuration->bottom_left_corner_e;

	// Rotation angle. Cos, sin, and tan are expensive computationally, so calculate once.
	rotation_angle = atan(east_width_m / north_height_m);

	cos_rotation_angle = cos(rotation_angle);
	sin_rotation_angle = sin(rotation_angle);

	total_height_m = sqrt(pow(cca_configuration->top_left_corner_n - cca_configuration->bottom_left_corner_n, 2.0f) +
						  pow(cca_configuration->top_left_corner_e - cca_configuration->bottom_left_corner_e, 2.0f));
	total_width_m  = sqrt(pow(cca_configuration->top_right_corner_n - cca_configuration->top_left_corner_n, 2.0f) +
						  pow(cca_configuration->top_right_corner_e - cca_configuration->top_left_corner_e, 2.0f));

        // Get the cos and sin for the Vs30 map rotation.
        cos_vs30_rotation_angle = cos(cca_vs30_map->rotation * DEG_TO_RAD);
        sin_vs30_rotation_angle = sin(cca_vs30_map->rotation * DEG_TO_RAD);

	// Let everyone know that we are initialized and ready for business.
	cca_is_initialized = 1;

	return SUCCESS;
}

/**
 * Queries CCA at the given points and returns the data that it finds.
 *
 * @param points The points at which the queries will be made.
 * @param data The data that will be returned (Vp, Vs, density, Qs, and/or Qp).
 * @param numpoints The total number of points to query.
 * @return SUCCESS or FAIL.
 */
int cca_query(cca_point_t *points, cca_properties_t *data, int numpoints) {
	int i = 0;
	double point_utm_e = 0, point_utm_n = 0;
	double temp_utm_e = 0, temp_utm_n = 0;
	int load_x_coord = 0, load_y_coord = 0, load_z_coord = 0;
	double x_percent = 0, y_percent = 0, z_percent = 0;
	cca_properties_t surrounding_points[8];

	int zone = 10;
	int longlat2utm = 0;

	for (i = 0; i < numpoints; i++) {

		// We need to be below the surface to service this query.
		if (points[i].depth < 0) {
			data[i].vp = -1;
			data[i].vs = -1;
			data[i].rho = -1;
			data[i].qp = -1;
			data[i].qs = -1;
			continue;
		}

		temp_utm_e = points[i].longitude; // * DEG_TO_RAD;
		temp_utm_n = points[i].latitude; // * DEG_TO_RAD;

		// Still need to use utm_geo
		utm_geo_(&temp_utm_e, &temp_utm_n, &point_utm_e, &point_utm_n, &zone, &longlat2utm);

		// Point within rectangle.
		point_utm_n -= cca_configuration->bottom_left_corner_n;
		point_utm_e -= cca_configuration->bottom_left_corner_e;

		temp_utm_n = point_utm_n;
		temp_utm_e = point_utm_e;

		// We need to rotate that point, the number of degrees we calculated above.
		point_utm_e = cos_rotation_angle * temp_utm_e - sin_rotation_angle * temp_utm_n;
		point_utm_n = sin_rotation_angle * temp_utm_e + cos_rotation_angle * temp_utm_n;

		// Which point base point does that correspond to?
		load_x_coord = floor(point_utm_e / total_width_m * (cca_configuration->nx - 1));
		load_y_coord = floor(point_utm_n / total_height_m * (cca_configuration->ny - 1));

		// And on the Z-axis?
		load_z_coord = (cca_configuration->depth / cca_configuration->depth_interval - 1) -
					   floor(points[i].depth / cca_configuration->depth_interval);

		// Are we outside the model's X and Y boundaries?
		if (load_x_coord > cca_configuration->nx - 2 || load_y_coord > cca_configuration->ny - 2 || load_x_coord < 0 || load_y_coord < 0) {
			data[i].vp = -1;
			data[i].vs = -1;
			data[i].rho = -1;
			data[i].qp = -1;
			data[i].qs = -1;
			continue;
		}

		// Get the X, Y, and Z percentages for the bilinear or trilinear interpolation below.
		x_percent = fmod(point_utm_e, total_width_m / cca_configuration->nx) / (total_width_m / cca_configuration->nx);
		y_percent = fmod(point_utm_n, total_height_m / cca_configuration->ny) / (total_height_m / cca_configuration->ny);
		z_percent = fmod(points[i].depth, cca_configuration->depth_interval) / cca_configuration->depth_interval;

		if (load_z_coord < 1) {
			// We're below the model boundaries. Bilinearly interpolate the bottom plane and use that value.
			data[i].vp = -1;
			data[i].vs = -1;
			data[i].rho = -1;
			data[i].qp = -1;
			data[i].qs = -1;

		        continue;
                } else {
if ((points[i].depth < cca_configuration->depth_interval) && (cca_configuration->gtl == 1)) {
                           get_vs30_based_gtl(&(points[i]), &(data[i]));
                           data[i].rho=calculate_density(data[i].vs);

                        } else {
			    // Read all the surrounding point properties.
			    read_properties(load_x_coord,     load_y_coord,     load_z_coord,     &(surrounding_points[0]));	// Orgin.
			    read_properties(load_x_coord + 1, load_y_coord,     load_z_coord,     &(surrounding_points[1]));	// Orgin + 1x
			    read_properties(load_x_coord,     load_y_coord + 1, load_z_coord,     &(surrounding_points[2]));	// Orgin + 1y
			    read_properties(load_x_coord + 1, load_y_coord + 1, load_z_coord,     &(surrounding_points[3]));	// Orgin + x + y, forms top plane.
			    read_properties(load_x_coord,     load_y_coord,     load_z_coord - 1, &(surrounding_points[4]));	// Bottom plane origin
			    read_properties(load_x_coord + 1, load_y_coord,     load_z_coord - 1, &(surrounding_points[5]));	// +1x
			    read_properties(load_x_coord,     load_y_coord + 1, load_z_coord - 1, &(surrounding_points[6]));	// +1y
			    read_properties(load_x_coord + 1, load_y_coord + 1, load_z_coord - 1, &(surrounding_points[7]));	// +x +y, forms bottom plane.

			    trilinear_interpolation(x_percent, y_percent, z_percent, surrounding_points, &(data[i]));
		}
}

		// Calculate Qp and Qs.
		if (data[i].vs < 1500)
			data[i].qs = data[i].vs * 0.02;
		else
			data[i].qs = data[i].vs * 0.10;

		data[i].qp = data[i].qs * 1.5;
	}

	return SUCCESS;
}

/**
 * Retrieves the material properties (whatever is available) for the given data point, expressed
 * in x, y, and z co-ordinates.
 *
 * @param x The x coordinate of the data point.
 * @param y The y coordinate of the data point.
 * @param z The z coordinate of the data point.
 * @param data The properties struct to which the material properties will be written.
 */
void read_properties(int x, int y, int z, cca_properties_t *data) {
	// Set everything to -1 to indicate not found.
	data->vp = -1;
	data->vs = -1;
	data->rho = -1;
	data->qp = -1;
	data->qs = -1;
	float *ptr = NULL;
	FILE *fp = NULL;
	int location = z * cca_configuration->nx * cca_configuration->ny + y * cca_configuration->nx + x;

	// Check our loaded components of the model.
	if (cca_velocity_model->vs_status == 2) {
		// Read from memory.
		ptr = (float *)cca_velocity_model->vs;
		data->vs = ptr[location];
	} else if (cca_velocity_model->vs_status == 1) {
		// Read from file.
		fp = (FILE *)cca_velocity_model->vs;
		fseek(fp, location * sizeof(float), SEEK_SET);
		fread(&(data->vs), sizeof(float), 1, fp);
	}

	// Check our loaded components of the model.
	if (cca_velocity_model->vp_status == 2) {
		// Read from memory.
		ptr = (float *)cca_velocity_model->vp;
		data->vp = ptr[location];
	} else if (cca_velocity_model->vp_status == 1) {
		// Read from file.
		fseek(fp, location * sizeof(float), SEEK_SET);
		fread(&(data->vp), sizeof(float), 1, fp);
	}

	// Check our loaded components of the model.
	if (cca_velocity_model->rho_status == 2) {
		// Read from memory.
		ptr = (float *)cca_velocity_model->rho;
		data->rho = ptr[location];
	} else if (cca_velocity_model->rho_status == 1) {
		// Read from file.
		fseek(fp, location * sizeof(float), SEEK_SET);
		fread(&(data->rho), sizeof(float), 1, fp);
	}
}

/**
 * Trilinearly interpolates given a x percentage, y percentage, z percentage and a cube of
 * data properties in top origin format (top plane first, bottom plane second).
 *
 * @param x_percent X percentage
 * @param y_percent Y percentage
 * @param z_percent Z percentage
 * @param eight_points Eight surrounding data properties
 * @param ret_properties Returned data properties
 */
void trilinear_interpolation(double x_percent, double y_percent, double z_percent,
							 cca_properties_t *eight_points, cca_properties_t *ret_properties) {
	cca_properties_t *temp_array = calloc(2, sizeof(cca_properties_t));
	cca_properties_t *four_points = eight_points;

	bilinear_interpolation(x_percent, y_percent, four_points, &temp_array[0]);

	// Now advance the pointer four "cca_properties_t" spaces.
	four_points += 4;

	// Another interpolation.
	bilinear_interpolation(x_percent, y_percent, four_points, &temp_array[1]);

	// Now linearly interpolate between the two.
	linear_interpolation(z_percent, &temp_array[0], &temp_array[1], ret_properties);

	free(temp_array);
}

/**
 * Bilinearly interpolates given a x percentage, y percentage, and a plane of data properties in
 * origin, bottom-right, top-left, top-right format.
 *
 * @param x_percent X percentage.
 * @param y_percent Y percentage.
 * @param four_points Data property plane.
 * @param ret_properties Returned data properties.
 */
void bilinear_interpolation(double x_percent, double y_percent, cca_properties_t *four_points, cca_properties_t *ret_properties) {
	cca_properties_t *temp_array = calloc(2, sizeof(cca_properties_t));
	linear_interpolation(x_percent, &four_points[0], &four_points[1], &temp_array[0]);
	linear_interpolation(x_percent, &four_points[2], &four_points[3], &temp_array[1]);
	linear_interpolation(y_percent, &temp_array[0], &temp_array[1], ret_properties);
	free(temp_array);
}

/**
 * Linearly interpolates given a percentage from x0 to x1, a data point at x0, and a data point at x1.
 *
 * @param percent Percent of the way from x0 to x1 (from 0 to 1 interval).
 * @param x0 Data point at x0.
 * @param x1 Data point at x1.
 * @param ret_properties Resulting data properties.
 */
void linear_interpolation(double percent, cca_properties_t *x0, cca_properties_t *x1, cca_properties_t *ret_properties) {
	ret_properties->vp  = (1 - percent) * x0->vp  + percent * x1->vp;
	ret_properties->vs  = (1 - percent) * x0->vs  + percent * x1->vs;
	ret_properties->rho = (1 - percent) * x0->rho + percent * x1->rho;
	ret_properties->qp  = (1 - percent) * x0->qp  + percent * x1->qp;
	ret_properties->qs  = (1 - percent) * x0->qs  + percent * x1->qs;
}

/**
 * Called when the model is being discarded. Free all variables.
 *
 * @return SUCCESS
 */
int cca_finalize() {
	pj_free(cca_latlon);
	pj_free(cca_utm);

	if (cca_velocity_model) free(cca_velocity_model);
	if (cca_configuration) free(cca_configuration);
        if (cca_vs30_map) free(cca_vs30_map);

	return SUCCESS;
}

/**
 * Returns the version information.
 *
 * @param ver Version string to return.
 * @param len Maximum length of buffer.
 * @return Zero
 */
int cca_version(char *ver, int len)
{
  int verlen;
  verlen = strlen(version_string);
  if (verlen > len - 1) {
    verlen = len - 1;
  }
  memset(ver, 0, len);
  strncpy(ver, version_string, verlen);
  return 0;
}

/**
 * Reads the cca_configuration file describing the various properties of CVM-S5 and populates
 * the cca_configuration struct. This assumes cca_configuration has been "calloc'ed" and validates
 * that each value is not zero at the end.
 *
 * @param file The cca_configuration file location on disk to read.
 * @param config The cca_configuration struct to which the data should be written.
 * @return Success or failure, depending on if file was read successfully.
 */
int read_cca_configuration(char *file, cca_configuration_t *config) {
	FILE *fp = fopen(file, "r");
	char key[40];
	char value[80];
	char line_holder[128];

	// If our file pointer is null, an error has occurred. Return fail.
	if (fp == NULL) {
		print_error("Could not open the cca_configuration file.");
		return FAIL;
	}

	// Read the lines in the cca_configuration file.
	while (fgets(line_holder, sizeof(line_holder), fp) != NULL) {
		if (line_holder[0] != '#' && line_holder[0] != ' ' && line_holder[0] != '\n') {
			sscanf(line_holder, "%s = %s", key, value);

			// Which variable are we editing?
			if (strcmp(key, "utm_zone") == 0) 				config->utm_zone = atoi(value);
			if (strcmp(key, "model_dir") == 0)				sprintf(config->model_dir, "%s", value);
			if (strcmp(key, "nx") == 0) 	  				config->nx = atoi(value);
			if (strcmp(key, "ny") == 0) 	  			 	config->ny = atoi(value);
			if (strcmp(key, "nz") == 0) 	  			 	config->nz = atoi(value);
			if (strcmp(key, "depth") == 0) 	  			 	config->depth = atof(value);
			if (strcmp(key, "top_left_corner_e") == 0) 		config->top_left_corner_e = atof(value);
			if (strcmp(key, "top_left_corner_n") == 0)		config->top_left_corner_n = atof(value);
			if (strcmp(key, "top_right_corner_e") == 0)		config->top_right_corner_e = atof(value);
			if (strcmp(key, "top_right_corner_n") == 0)		config->top_right_corner_n = atof(value);
			if (strcmp(key, "bottom_left_corner_e") == 0)	config->bottom_left_corner_e = atof(value);
			if (strcmp(key, "bottom_left_corner_n") == 0)	config->bottom_left_corner_n = atof(value);
			if (strcmp(key, "bottom_right_corner_e") == 0)	config->bottom_right_corner_e = atof(value);
			if (strcmp(key, "bottom_right_corner_n") == 0)	config->bottom_right_corner_n = atof(value);
			if (strcmp(key, "depth_interval") == 0)			config->depth_interval = atof(value);
                        if (strcmp(key, "p0") == 0)                                             config->p0 = atof(value);
                        if (strcmp(key, "p1") == 0)                                             config->p1 = atof(value);
                        if (strcmp(key, "p2") == 0)                                             config->p2 = atof(value);
                        if (strcmp(key, "p3") == 0)                                             config->p3 = atof(value);
                        if (strcmp(key, "p4") == 0)                                             config->p4 = atof(value);
                        if (strcmp(key, "p5") == 0)                                             config->p5 = atof(value);
                        if (strcmp(key, "gtl") == 0) {
                                if (strcmp(value, "on") == 0) config->gtl = 1;
                                else config->gtl = 0;
                        }
		}
	}

	// Have we set up all cca_configuration parameters?
	if (config->utm_zone == 0 || config->nx == 0 || config->ny == 0 || config->nz == 0 || config->model_dir[0] == '\0' ||
		config->top_left_corner_e == 0 || config->top_left_corner_n == 0 || config->top_right_corner_e == 0 ||
		config->top_right_corner_n == 0 || config->bottom_left_corner_e == 0 || config->bottom_left_corner_n == 0 ||
		config->bottom_right_corner_e == 0 || config->bottom_right_corner_n == 0 || config->depth == 0 ||
		config->depth_interval == 0) {
		print_error("One cca_configuration parameter not specified. Please check your cca_configuration file.");
		return FAIL;
	}

	fclose(fp);

	return SUCCESS;
}

/**
 *  * Calculates the density based off of Vs. Based on Nafe-Drake scaling relationship.
 *   *
 *    * @param vs The Vs value off which to scale.
 *     * @return Density, in g/m^3.
 *      */
double calculate_density(double vs) {
        double retVal;
        vs = vs / 1000;
        retVal = cca_configuration->p0 + cca_configuration->p1 * vs + cca_configuration->p2 * pow(vs, 2) +
                         cca_configuration->p3 * pow(vs, 3) + cca_configuration->p4 * pow(vs, 4) + cca_configuration->p5 * pow(vs, 5);
        retVal = retVal * 1000;
        return retVal;
}


/**
 * Prints the error string provided.
 *
 * @param err The error string to print out to stderr.
 */
void print_error(char *err) {
	fprintf(stderr, "An error has occurred while executing CCA. The error was:\n\n");
	fprintf(stderr, "%s", err);
	fprintf(stderr, "\n\nPlease contact software@scec.org and describe both the error and a bit\n");
	fprintf(stderr, "about the computer you are running CCA on (Linux, Mac, etc.).\n");
}

/**
 * Tries to read the model into memory.
 *
 * @param model The model parameter struct which will hold the pointers to the data either on disk or in memory.
 * @return 2 if all files are read to memory, SUCCESS if file is found but at least 1
 * is not in memory, FAIL if no file found.
 */
int try_reading_model(cca_model_t *model) {
	double base_malloc = cca_configuration->nx * cca_configuration->ny * cca_configuration->nz * sizeof(float);
	int file_count = 0;
	int all_read_to_memory = 1;
	char current_file[128];
	FILE *fp;

	// Let's see what data we actually have.
	sprintf(current_file, "%s/vp.dat", cca_iteration_directory);
	if (access(current_file, R_OK) == 0) {
		model->vp = malloc(base_malloc);
		if (model->vp != NULL) {
			// Read the model in.
			fp = fopen(current_file, "rb");
			fread(model->vp, 1, base_malloc, fp);
			fclose(fp);
			model->vp_status = 2;
		} else {
			all_read_to_memory = 0;
			model->vp = fopen(current_file, "rb");
			model->vp_status = 1;
		}
		file_count++;
	}

	sprintf(current_file, "%s/vs.dat", cca_iteration_directory);
	if (access(current_file, R_OK) == 0) {
		model->vs = malloc(base_malloc);
		if (model->vs != NULL) {
			// Read the model in.
			fp = fopen(current_file, "rb");
			fread(model->vs, 1, base_malloc, fp);
			fclose(fp);
			model->vs_status = 2;
		} else {
			all_read_to_memory = 0;
			model->vs = fopen(current_file, "rb");
			model->vs_status = 1;
		}
		file_count++;
	}

	sprintf(current_file, "%s/density.dat", cca_iteration_directory);
	if (access(current_file, R_OK) == 0) {
		model->rho = malloc(base_malloc);
		if (model->rho != NULL) {
			// Read the model in.
			fp = fopen(current_file, "rb");
			fread(model->rho, 1, base_malloc, fp);
			fclose(fp);
			model->rho_status = 2;
		} else {
			all_read_to_memory = 0;
			model->rho = fopen(current_file, "rb");
			model->rho_status = 1;
		}
		file_count++;
	}

	sprintf(current_file, "%s/qp.dat", cca_iteration_directory);
	if (access(current_file, R_OK) == 0) {
		model->qp = malloc(base_malloc);
		if (model->qp != NULL) {
			// Read the model in.
			fp = fopen(current_file, "rb");
			fread(model->qp, 1, base_malloc, fp);
			fclose(fp);
			model->qp_status = 2;
		} else {
			all_read_to_memory = 0;
			model->qp = fopen(current_file, "rb");
			model->qp_status = 1;
		}
		file_count++;
	}

	sprintf(current_file, "%s/qs.dat", cca_iteration_directory);
	if (access(current_file, R_OK) == 0) {
		model->qs = malloc(base_malloc);
		if (model->qs != NULL) {
			// Read the model in.
			fp = fopen(current_file, "rb");
			fread(model->qs, 1, base_malloc, fp);
			fclose(fp);
			model->qs_status = 2;
		} else {
			all_read_to_memory = 0;
			model->qs = fopen(current_file, "rb");
			model->qs_status = 1;
		}
		file_count++;
	}

	if (file_count == 0)
		return FAIL;
	else if (file_count > 0 && all_read_to_memory == 0)
		return SUCCESS;
	else
		return 2;
}


/**
 * Reads the format of the Vs30 data e-tree. This file location is typically specified
 * in the cca_configuration file of the model.
 *
 * @param filename The e-tree's file location from which to read.
 * @param map The outputted map cca_configuration structure.
 */
int read_vs30_map(char *filename, cca_vs30_map_config_t *map) {
	char appmeta[512];
	char *token;
	int index = 0, retVal = 0;
	map->vs30_map = etree_open(filename, O_RDONLY, 64, 0, 3);
	retVal = snprintf(appmeta, sizeof(appmeta), "%s", etree_getappmeta(map->vs30_map));

	if (retVal >= 0 && retVal < 128) {
		return FAIL;
	}

	// Now we need to parse the map cca_configuration.
	index = 0;
	token = strtok(appmeta, "|");

	while (token != NULL) {
		switch (index) {
	    case 0:
	    	snprintf(map->type, sizeof(map->type), "%s", token);
	    	break;
	    case 1:
	    	snprintf(map->description, sizeof(map->description), "%s", token);
	    	break;
	    case 2:
	    	snprintf(map->author, sizeof(map->author), "%s", token);
	    	break;
	    case 3:
	    	snprintf(map->date, sizeof(map->date), "%s", token);
	    	break;
	    case 4:
	    	sscanf(token, "%lf", &(map->spacing));
	    	break;
	    case 5:
	    	snprintf(map->schema, sizeof(map->schema), "%s", token);
	    	break;
	    case 6:
	    	snprintf(map->projection, sizeof(map->projection), "%s", token);
	    	break;
	    case 7:
	    	sscanf(token, "%lf,%lf,%lf", &(map->origin_point.longitude), &(map->origin_point.latitude),
	    			&(map->origin_point.depth));
	    	break;
	    case 8:
	    	sscanf(token, "%lf", &(map->rotation));
	    	break;
	    case 9:
	    	sscanf(token, "%lf,%lf,%lf", &(map->x_dimension), &(map->y_dimension), &(map->z_dimension));
	    	break;
	    case 10:
	    	sscanf(token, "%u,%u,%u", &(map->x_ticks), &(map->y_ticks), &(map->z_ticks));
	    	break;
	    default:
	    	fprintf(stderr, "Unexpected metadata. Please check your Vs30 e-tree within UCVM.\n");
	    	return FAIL;
	    	break;
		}
	    index++;
	    token = strtok(NULL, "|");
	}

	return SUCCESS;

}

/**
 * Given a latitude and longitude in WGS84 co-ordinates, we find the corresponding e-tree octant
 * in the Vs30 map e-tree and read the value as well as interpolate bilinearly.
 *
 * @param longitude The longitude in WGS84 format.
 * @param latitude The latitude in WGS84 format.
 * @param map The Vs30 map structure as defined during the initialization procedure.
 * @return The Vs30 value at that point, or -1 if outside the boundaries.
 */
double get_vs30_value(double longitude, double latitude, cca_vs30_map_config_t *map) {
	// Convert both points to UTM.
	double longitude_utm_e = longitude * DEG_TO_RAD;
	double latitude_utm_n = latitude * DEG_TO_RAD;
	double vs30_long_utm_e = map->origin_point.longitude * DEG_TO_RAD;
	double vs30_lat_utm_n = map->origin_point.latitude * DEG_TO_RAD;
	double temp_rotated_point_n = 0.0, temp_rotated_point_e = 0.0;
	double rotated_point_n = 0.0, rotated_point_e = 0.0;
	double percent = 0.0;
	int loc_x = 0, loc_y = 0;
	etree_addr_t addr;
	vs30_mpayload_t vs30_payload[4];

	int max_level = ceil(log(map->x_dimension / map->spacing) / log(2.0));

	etree_tick_t edgetics = (etree_tick_t)1 << (ETREE_MAXLEVEL - max_level);
	double map_edgesize = map->x_dimension / (double)((etree_tick_t)1<<max_level);

	pj_transform(cca_latlon, cca_aeqd, 1, 1, &longitude_utm_e, &latitude_utm_n, NULL);
	pj_transform(cca_latlon, cca_aeqd, 1, 1, &vs30_long_utm_e, &vs30_lat_utm_n, NULL);

	// Now that both are in UTM, we can subtract and rotate.
	temp_rotated_point_e = longitude_utm_e - vs30_long_utm_e;
	temp_rotated_point_n = latitude_utm_n - vs30_lat_utm_n;

	rotated_point_e = cos_vs30_rotation_angle * temp_rotated_point_e - sin_vs30_rotation_angle * temp_rotated_point_n;
	rotated_point_n = sin_vs30_rotation_angle * temp_rotated_point_e + cos_vs30_rotation_angle * temp_rotated_point_n;

	// Are we within the box?
	if (rotated_point_e < 0 || rotated_point_n < 0 || rotated_point_e > map->x_dimension ||
		rotated_point_n > map->y_dimension) return -1;

	// Get the integer location of the grid point within the map.
	loc_x = floor(rotated_point_e / map_edgesize);
	loc_y = floor(rotated_point_n / map_edgesize);

	// We need the four surrounding points for bilinear interpolation.
	addr.level = ETREE_MAXLEVEL;
	addr.x = loc_x * edgetics; addr.y = loc_y * edgetics; addr.z = 0;
    /* Adjust addresses for edges of grid */
    if (addr.x >= map->x_ticks) addr.x = map->x_ticks - edgetics;
    if (addr.y >= map->y_ticks) addr.y = map->y_ticks - edgetics;
	etree_search(map->vs30_map, addr, NULL, "*", &(vs30_payload[0]));
	addr.x = (loc_x + 1) * edgetics; addr.y = loc_y * edgetics;
    if (addr.x >= map->x_ticks) addr.x = map->x_ticks - edgetics;
    if (addr.y >= map->y_ticks) addr.y = map->y_ticks - edgetics;
	etree_search(map->vs30_map, addr, NULL, "*", &(vs30_payload[1]));
	addr.x = loc_x * edgetics; addr.y = (loc_y + 1) * edgetics;
    if (addr.x >= map->x_ticks) addr.x = map->x_ticks - edgetics;
    if (addr.y >= map->y_ticks) addr.y = map->y_ticks - edgetics;
	etree_search(map->vs30_map, addr, NULL, "*", &(vs30_payload[2]));
	addr.x = (loc_x + 1) * edgetics; addr.y = (loc_y + 1) * edgetics;
    if (addr.x >= map->x_ticks) addr.x = map->x_ticks - edgetics;
    if (addr.y >= map->y_ticks) addr.y = map->y_ticks - edgetics;
	etree_search(map->vs30_map, addr, NULL, "*", &(vs30_payload[3]));

	percent = fmod(rotated_point_e / map->spacing, map->spacing) / map->spacing;
	vs30_payload[0].vs30 = percent * vs30_payload[0].vs30 + (1 - percent) * vs30_payload[1].vs30;
	vs30_payload[1].vs30 = percent * vs30_payload[2].vs30 + (1 - percent) * vs30_payload[3].vs30;

	return vs30_payload[0].vs30;
}

/**
 * Gets the GTL value using the Wills and Wald dataset, given a latitude, longitude and depth.
 *
 * @param point The point at which to retrieve the property. Note, depth is ignored.
 * @param data The material properties at the point specified, or -1 if not found.
 * @return Success or failure.
 */
int get_vs30_based_gtl(cca_point_t *point, cca_properties_t *data) {
        double a = 0.5, b = 0.6, c = 0.5;
	double percent_z = point->depth / cca_configuration->depth_interval;
	double f = 0.0, g = 0.0;
	double vs30 = 0.0, vp30 = 0.0;

	// Double check that we're above the first layer.
	if (percent_z > 1) return FAIL;

	// Query for the point at depth_interval.
	cca_point_t *pt = calloc(1, sizeof(cca_point_t));
	cca_properties_t *dt = calloc(1, sizeof(cca_properties_t));

	pt->latitude = point->latitude;
	pt->longitude = point->longitude;
	pt->depth = cca_configuration->depth_interval;

	if (cca_query(pt, dt, 1) != SUCCESS) return FAIL;

	// Now we need the Vs30 data value.
	vs30 = get_vs30_value(point->longitude, point->latitude, cca_vs30_map);

	if (vs30 == -1) {
		data->vp = -1;
		data->vs = -1;
	} else {
		// Get the point's material properties within the GTL.
		f = percent_z + b * (percent_z - pow(percent_z, 2.0f));
		g = a - a * percent_z + c * (pow(percent_z, 2.0f) + 2.0 * sqrt(percent_z) - 3.0 * percent_z);
		data->vs = f * dt->vs + g * vs30;
//fprintf(stderr,"XXX f %f and g %f\n", f, g);
		vs30 = vs30 / 1000;
		vp30 = 0.9409 + 2.0947 * vs30 - 0.8206 * pow(vs30, 2.0f) + 0.2683 * pow(vs30, 3.0f) - 0.0251 * pow(vs30, 4.0f);
		vp30 = vp30 * 1000;
		data->vp = f * dt->vp + g * vp30;
	}

	free(pt);
	free(dt);

	return SUCCESS;
}


// The following functions are for dynamic library mode. If we are compiling
// a static library, these functions must be disabled to avoid conflicts.
#ifdef DYNAMIC_LIBRARY

/**
 * Init function loaded and called by the UCVM library. Calls cca_init.
 *
 * @param dir The directory in which UCVM is installed.
 * @return Success or failure.
 */
int model_init(const char *dir, const char *label) {
	return cca_init(dir, label);
}

/**
 * Query function loaded and called by the UCVM library. Calls cca_query.
 *
 * @param points The basic_point_t array containing the points.
 * @param data The basic_properties_t array containing the material properties returned.
 * @param numpoints The number of points in the array.
 * @return Success or fail.
 */
int model_query(cca_point_t *points, cca_properties_t *data, int numpoints) {
	return cca_query(points, data, numpoints);
}

/**
 * Finalize function loaded and called by the UCVM library. Calls cca_finalize.
 *
 * @return Success
 */
int model_finalize() {
	return cca_finalize();
}

/**
 * Version function loaded and called by the UCVM library. Calls cca_version.
 *
 * @param ver Version string to return.
 * @param len Maximum length of buffer.
 * @return Zero
 */
int model_version(char *ver, int len) {
	return cca_version(ver, len);
}

int (*get_model_init())(const char *, const char *) {
        return &cca_init;
}
int (*get_model_query())(cca_point_t *, cca_properties_t *, int) {
         return &cca_query;
}
int (*get_model_finalize())() {
         return &cca_finalize;
}
int (*get_model_version())(char *, int) {
         return &cca_version;
}

#endif
