//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//	SPDX-License-Identifier: BSD-3-Clause

// =============================================================================
// This demo presents a cone penetrameter test with a soil sample made of clumped
// particles of various sizes. Before the test starts, when compress the terrain
// first, and note that the compressor used in this process has its position
// explicitly controlled step-by-step.
// =============================================================================

#include <DEM/API.h>
#include <DEM/HostSideHelpers.hpp>
#include <DEM/utils/Samplers.hpp>
#include <Utils.h>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <map>
#include <random>

using namespace deme;

const double math_PI = 3.14159;

enum FAMILY : int {
    FREE = 0,
    TRANSLATING,
    FIXED
};

//Function that uses inputted data_dir and file label to create new directory for results
std::filesystem::path InitializeOutputDirectories(std::filesystem::path data_dir, std::string label) {
    std::filesystem::path out_dir = data_dir / label / Utils::getCurrentTimeStamp();

    std::error_code ec;
    if (!std::filesystem::create_directories(out_dir, ec)) {
        throw std::runtime_error("Failed to create output directories.");
    }

    std::cout << "Output directory: " << out_dir << std::endl;

    return out_dir;
}

int main(int argc, char* argv[]) {

    if (argc != 5) {
        std::cerr << "Usage: ./SoilPenetrometer <scale_factor> <terrain_file_path> <data_dir> <label>" << std::endl;
        return EXIT_FAILURE;
    }

    double scale_factor = std::atof(argv[1]);
    std::filesystem::path terrain_filepath = argv[2];
    std::filesystem::path data_dir = argv[3];
    std::string label = argv[4];

    // Calls function to create output directory
    std::filesystem::path out_dir = InitializeOutputDirectories(data_dir, label);

    DEMSolver DEMSim;
    DEMSim.SetVerbosity(INFO);
    DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV);
    DEMSim.SetOutputContent(OUTPUT_CONTENT::ABSV);
    DEMSim.SetMeshOutputFormat(MESH_FORMAT::VTK);
    DEMSim.SetContactOutputContent({"OWNER", "FORCE", "POINT"});

    // E, nu, CoR, mu, Crr...
    auto mat_type_cone = DEMSim.LoadMaterial({{"E", 1e9}, {"nu", 0.3}, {"CoR", 0.8}, {"mu", 0.7}, {"Crr", 0.00}});
    auto mat_type_terrain = DEMSim.LoadMaterial({{"E", 1e9}, {"nu", 0.3}, {"CoR", 0.8}, {"mu", 0.4}, {"Crr", 0.00}});
    // If you don't have this line, then values will take average between 2 materials, when they are in contact
    DEMSim.SetMaterialPropertyPair("CoR", mat_type_cone, mat_type_terrain, 0.8);
    DEMSim.SetMaterialPropertyPair("mu", mat_type_cone, mat_type_terrain, 0.7);

    float cone_speed = 0.03;
    float step_size = 5e-6;
    double soil_bin_diameter = 0.584;
    double cone_surf_area = 323e-6;
    double cone_diameter = std::sqrt(cone_surf_area / math_PI) * 2;
    double bottom = -0.5;
    

    //world size and particle initialization
    {

    std::cout << "Defining World..." << std::endl;

    // Define world dimensions
    double world_size_x = 1.0;
    double world_size_y = 0.3;
    double world_size_z = 2.0;
    DEMSim.InstructBoxDomainDimension(world_size_x, world_size_y, world_size_z);
    DEMSim.InstructBoxDomainBoundingBC("top_open", mat_type_terrain);

    // Add boundary planes
    float bottom = -0.5f;
    DEMSim.AddBCPlane(make_float3(0.0f, 0.0f, bottom), make_float3(0.0f, 0.0f, 1.0f), mat_type_terrain);
    DEMSim.AddBCPlane(make_float3(0.0f, static_cast<float>(world_size_y) / 2.0f, 0.0f), make_float3(0.0f, -1.0f, 0.0f), mat_type_terrain);
    DEMSim.AddBCPlane(make_float3(0.0f, -static_cast<float>(world_size_y) / 2.0f, 0.0f), make_float3(0.0f, 1.0f, 0.0f), mat_type_terrain);
    DEMSim.AddBCPlane(make_float3(-static_cast<float>(world_size_x) / 2.0f, 0.0f, 0.0f), make_float3(1.0f, 0.0f, 0.0f), mat_type_terrain);
    DEMSim.AddBCPlane(make_float3(static_cast<float>(world_size_x) / 2.0f, 0.0f, 0.0f), make_float3(-1.0f, 0.0f, 0.0f), mat_type_terrain);

    // Define terrain particle templates
    float terrain_density = 2.6e3f;
    float volume1 = 4.2520508f;
    float mass1 = terrain_density * volume1;
    float3 MOI1 = make_float3(1.6850426f, 1.6375114f, 2.1187753f) * terrain_density;
    float volume2 = 2.1670011f;
    float mass2 = terrain_density * volume2;
    float3 MOI2 = make_float3(0.57402126f, 0.60616378f, 0.92890173f) * terrain_density;


    // Scale factors
    std::vector<double> scales = {0.0014, 0.00075833, 0.00044, 0.0003, 0.0002, 0.00018333, 0.00017};
    for (auto& scale : scales) {
        scale *= scale_factor * 10.0;
    }

    std::cout << "Loading clump templates..." << std::endl;


    // Load clump templates
    std::shared_ptr<DEMClumpTemplate> my_template2 = DEMSim.LoadClumpType(mass2, MOI2, GetDEMEDataFile("clumps/triangular_flat_6comp.csv"), mat_type_terrain);
    std::shared_ptr<DEMClumpTemplate> my_template1 = DEMSim.LoadClumpType(mass1, MOI1, GetDEMEDataFile("clumps/triangular_flat.csv"), mat_type_terrain);
    std::vector<std::shared_ptr<DEMClumpTemplate>> ground_particle_templates = {
        my_template2,
        DEMSim.Duplicate(my_template2),
        my_template1,
        DEMSim.Duplicate(my_template1),
        DEMSim.Duplicate(my_template1),
        DEMSim.Duplicate(my_template1),
        DEMSim.Duplicate(my_template1)
    };

    // Scale and name templates
    for (size_t i = 0; i < scales.size(); ++i) {
        auto& tmpl = ground_particle_templates.at(i);
        tmpl->Scale(scales.at(i));

        char t_name[20];
        std::sprintf(t_name, "%04zu", i);
        tmpl->AssignName(std::string(t_name));
    }

    // Load clump locations from file
    std::cout << "Making terrain..." << std::endl;
    std::unordered_map<std::string, std::vector<float3>> clump_xyz;
    std::unordered_map<std::string, std::vector<float4>> clump_quaternion;

    try {
        clump_xyz = DEMSim.ReadClumpXyzFromCsv(terrain_filepath);
        clump_quaternion = DEMSim.ReadClumpQuatFromCsv(terrain_filepath);
    } catch (...) {
        throw std::runtime_error("Failed to read clump checkpoint file. Ensure the file exists and is correctly formatted.");
    }

    std::vector<float3> in_xyz;
    std::vector<float4> in_quat;
    std::vector<std::shared_ptr<DEMClumpTemplate>> in_types;
    unsigned int t_num = 0;

    for (const auto& scale : scales) {
        char t_name[20];
        std::sprintf(t_name, "%04u", t_num);

        auto this_type_xyz = clump_xyz[std::string(t_name)];
        auto this_type_quat = clump_quaternion[std::string(t_name)];

        size_t n_clump_this_type = this_type_xyz.size();
        std::cout << "Loading clump " << std::string(t_name) << " with " << n_clump_this_type << " particles." << std::endl;

        std::vector<std::shared_ptr<DEMClumpTemplate>> this_type(n_clump_this_type, ground_particle_templates.at(t_num));

        in_xyz.insert(in_xyz.end(), this_type_xyz.begin(), this_type_xyz.end());
        in_quat.insert(in_quat.end(), this_type_quat.begin(), this_type_quat.end());
        in_types.insert(in_types.end(), this_type.begin(), this_type.end());

        std::cout << "Added clump type " << t_num << std::endl;
        t_num++;
    }

    // Create and add clump batch
    DEMClumpBatch base_batch(in_xyz.size());
    base_batch.SetTypes(in_types);
    base_batch.SetPos(in_xyz);
    base_batch.SetOriQ(in_quat);

    DEMSim.AddClumps(base_batch);
    }


    // Load in the cone used for this penetration test
    auto cone_tip = DEMSim.AddWavefrontMeshObject(GetDEMEDataFile("mesh/cone.obj"), mat_type_cone);
    float tip_height = std::sqrt(3.);
    std::shared_ptr<DEMTracker> tip_tracker;
    std::shared_ptr<DEMTracker> body_tracker;
    {
        auto cone_body = DEMSim.AddWavefrontMeshObject(GetDEMEDataFile("mesh/cyl_r1_h2.obj"), mat_type_cone);
        std::cout << "Total num of triangles: " << cone_tip->GetNumTriangles() + cone_body->GetNumTriangles() << std::endl;

        // The initial cone mesh has base radius 1, and height 1. Let's stretch it a bit so it has a 60deg tip, instead of
        // 90deg.
        cone_tip->Scale(make_float3(1, 1, tip_height));
        // Then set mass properties
        float cone_mass = 7.8e3 * tip_height / 3 * math_PI;
        cone_tip->SetMass(cone_mass);
        // You can checkout https://en.wikipedia.org/wiki/List_of_moments_of_inertia
        cone_tip->SetMOI(make_float3(cone_mass * (3. / 20. + 3. / 80. * tip_height * tip_height),
                                    cone_mass * (3. / 20. + 3. / 80. * tip_height * tip_height), 3 * cone_mass / 10));
        // This cone mesh has its tip at the origin. And, float4 quaternion pattern is (x, y, z, w).
        cone_tip->InformCentroidPrincipal(make_float3(0, 0, 3. / 4. * tip_height), make_float4(0, 0, 0, 1));
        // Note the scale method will scale mass and MOI automatically. But this only goes for the case you scale xyz all
        // together; otherwise, the MOI scaling will not be accurate and you should manually reset them.
        cone_tip->Scale(cone_diameter / 2);
        cone_tip->SetFamily(FAMILY::FIXED);

        // The define the body that is connected to the tip
        float body_mass = 7.8e3 * math_PI;
        cone_body->SetMass(body_mass);
        cone_body->SetMOI(make_float3(body_mass * 7 / 12, body_mass * 7 / 12, body_mass / 2));
        // This cyl mesh (h = 2m, r = 1m) has its center at the origin. So the following call actually has no effect...
        cone_body->InformCentroidPrincipal(make_float3(0, 0, 0), make_float4(0, 0, 0, 1));
        cone_body->Scale(make_float3(cone_diameter / 2, cone_diameter / 2, 0.5));
        cone_body->SetFamily(FAMILY::FIXED);

        // Track the cone_tip
        tip_tracker = DEMSim.Track(cone_tip);
        body_tracker = DEMSim.Track(cone_body);
    }

    // Because the cone's motion is completely pre-determined, we can just prescribe family 1
    DEMSim.SetFamilyPrescribedLinVel(FAMILY::TRANSLATING, "0", "0", "-" + to_string_with_precision(cone_speed));
    // Cone is initially in family 2, sleeping...
    DEMSim.SetFamilyFixed(FAMILY::FIXED);
    DEMSim.DisableContactBetweenFamilies(FAMILY::FREE, FAMILY::FIXED);

    // Some inspectors
    auto max_z_finder = DEMSim.CreateInspector("clump_max_z");
    auto total_mass_finder = DEMSim.CreateInspector("clump_mass");

    //Initializing DEM Sim
    {
        DEMSim.SetInitTimeStep(step_size);
        DEMSim.SetGravitationalAcceleration(make_float3(0, 0, -9.81));
        // CD freq will be auto-adapted so it does not matter much here.
        DEMSim.SetCDUpdateFreq(20);
        // Max velocity info is generally just for the solver's reference and the user do not have to set it. The solver
        // wouldn't take into account a vel larger than this when doing async-ed contact detection: but this vel won't
        // happen anyway and if it does, something already went wrong.
        DEMSim.SetMaxVelocity(10.);
        DEMSim.SetFamilyFixed(10);
        DEMSim.Initialize();
    }

    unsigned int fps = 20;
    float terrain_max_z = max_z_finder->GetValue();
    float bulk_density;
    unsigned int currframe = 0;
    
    // Code to set precision of scale_factor when converting to a string
    std::ostringstream stream;
    stream << std::setprecision(5) << scale_factor;
    std::string scale = stream.str();

    std::filesystem::path output_datafile_path = out_dir / "output.csv";
    std::cout << output_datafile_path << std::endl;
    std::ofstream output_datafile;
    output_datafile.open(output_datafile_path);
    if (!output_datafile) {
        throw std::runtime_error("Failed to open output.csv for writing.");
    }
    output_datafile << "t,pos_z,pressure" <<std::endl;

    float sim_end = 7.0;
    fps = 2500;
    float frame_time = 1.0 / fps;
    std::cout << "Output at " << fps << " FPS" << std::endl;

    // Put the cone in place
    double starting_height = terrain_max_z + 0.03;
    // Its initial position should be right above the cone tip...
    body_tracker->SetPos(make_float3(0, 0, 0.5 + (cone_diameter / 2 / 4 * tip_height) + starting_height));
    // Note that position of objects is always the location of their centroid
    tip_tracker->SetPos(make_float3(0, 0, starting_height));
    // The tip location, used to measure penetration length
    double tip_starting_z = -cone_diameter / 2 * 3 / 4 * tip_height + starting_height;
    double tip_z = tip_starting_z;

    // Enable cone
    DEMSim.ChangeFamily(FAMILY::FIXED, FAMILY::TRANSLATING);

    double tip_z_when_first_hit;
    bool hit_terrain = false;
    unsigned int frame_count = 0;

    float box_halfsize_x = 0.125;
    float box_halfsize_y = 0.125;

    DEMSim.DoDynamicsThenSync(0.0f);
    DEMSim.ChangeClumpFamily(FAMILY::FIXED);

    float pos_x = 0;
    float pos_y = 0;

    std::pair<float, float> Xrange = {pos_x - box_halfsize_x, pos_x + box_halfsize_x};
    std::pair<float, float> Yrange = {pos_y - box_halfsize_y, pos_y + box_halfsize_y};
    DEMSim.ChangeClumpFamily(FAMILY::FREE, Xrange, Yrange);
    
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    for (float t = 0; t < sim_end; t += frame_time) {

        float3 forces = tip_tracker->ContactAcc();
        // Note cone_mass is not the true mass, b/c we scaled the the cone tip! So we use true mass by using
        // cone_tip->mass.
        forces *= cone_tip->mass;
        float pressure = std::abs(forces.z) / cone_surf_area;
        if (pressure > 1e-4 && !hit_terrain) {
            hit_terrain = true;
            tip_z_when_first_hit = tip_z;
        }
        float penetration = (hit_terrain) ? tip_z_when_first_hit - tip_z : 0;
        std::cout << "Time: " << t << std::endl;
        std::cout << "Z coord of tip: " << tip_z << std::endl;
        std::cout << "Penetration: " << penetration << std::endl;
        std::cout << "Force on cone: " << forces.x << ", " << forces.y << ", " << forces.z << std::endl;
        std::cout << "Pressure: " << pressure << std::endl;

        //Writing to output.csv file
        try {
            output_datafile<< t << ","
                       << tip_z << ","
                       << pressure
                       << std::endl;
            output_datafile.flush();
        } catch (...) {
            throw std::runtime_error("Unable to write content to output file.");
        }

        float bottom_buffer = 0.01;
        if (frame_count % 500 == 0) {
            char filename[200], meshname[200];
            std::cout << "Outputting frame: " << currframe << std::endl;
            sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), currframe);
            sprintf(meshname, "%s/DEMdemo_mesh_%04d.vtk", out_dir.c_str(), currframe++);
            DEMSim.WriteSphereFile(std::string(filename));
            DEMSim.WriteMeshFile(std::string(meshname));
            DEMSim.ShowThreadCollaborationStats();
        }

        //Termination Condition
        if ((tip_z) < (bottom + bottom_buffer)) {
        std::cout << "This is far enough, stopping the simulation..." << std::endl;
        DEMSim.DoDynamicsThenSync(0.0f);
        break;
        }
        
        DEMSim.DoDynamicsThenSync(frame_time);
        tip_z -= cone_speed * frame_time;

        frame_count++;
    }
    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_sec = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    std::cout << time_sec.count() << " seconds (wall time) to finish the simulation" << std::endl;

    std::cout << "ConePenetration demo exiting..." << std::endl;
    return 0;
}