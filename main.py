def main():
    from Ores import Ore
    from CourseGrinder import CoarseCrusher
    from Washer import Washer
    from CrusherGrinder import CrusherGrinder
    from Dryer import Dryer
    from Calcinator import Calcinator
    from PretreatmentPipeline import PretreatmentPipeline

    print("=== HAIMAN Process Pretreatment (Continuous) ===")

    # 1) Ask user: do we want to specify total mass or total run time?
    mode_input = input("Type '1' to specify total mass, '2' for run time: ")
    if mode_input == "2":
        # user sets operation_time (hours) and throughput (kg/h)
        time_input = input("Enter total operation time (hours): ")
        try:
            operation_time = float(time_input)
        except ValueError:
            operation_time = 10.0

        feed_input = input("Enter throughput (kg/h): ")
        try:
            throughput = float(feed_input)
        except ValueError:
            throughput = 1000.0

        total_ore_mass = operation_time * throughput
        print(f"\nWe will process {total_ore_mass:.2f} kg of ore over {operation_time} h at {throughput} kg/h.\n")

    else:
        # user sets total mass, plus a throughput => we derive total run time
        mass_input = input("Enter total mass of ore (kg): ")
        try:
            total_ore_mass = float(mass_input)
        except ValueError:
            total_ore_mass = 1000.0

        feed_input = input("Enter throughput (kg/h) [default=1000]: ")
        try:
            throughput = float(feed_input)
        except ValueError:
            throughput = 1000.0

        # We'll compute final time after we run the pipeline
        print(f"\nWe have {total_ore_mass:.2f} kg of ore, throughput {throughput} kg/h.\n")

    # 2) Get initial temperature, density, etc. (like before)
    temp_input = input("Enter initial ore temperature (°C) [default 25]: ")
    try:
        ore_temp = float(temp_input)
    except ValueError:
        ore_temp = 25.0

    dens_input = input("Enter initial bulk density (kg/m^3) [default 2800]: ")
    try:
        ore_density = float(dens_input)
    except ValueError:
        ore_density = 2800.0

    # 3) Define phases (sum=1.0) and size distribution
    phases = {
        "MnO": 0.1657,
        "MnO2": 0.1655,
        "Mn2O3": 0.1908,
        "Mn3O4": 0.0,
        "Fe2O3": 0.13,
        "SiO2": 0.052,
        "Al2O3": 0.0017,
        "CaO": 0.121,
        "MgO": 0.032,
        "C": 0.001,
        "H2O": 0.0015,
        "Impurities": 0.1388,
    }

    size_dist = {
        ">50 mm": 0.80,
        "9.5-50 mm": 0.1,
        "2.5-9.5 mm": 0.1,
        "<2.5 mm": 0
    }

    ore = Ore(
        mass=total_ore_mass,
        temperature=ore_temp,
        density=ore_density,
        phases=phases,
        size_distribution=size_dist
    )

    # 4) Phase densities (optional)
    PHASE_DENSITIES = {
        "MnO": 5037,
        "MnO2": 5030,
        "Mn2O3": 4380,
        "Mn3O4": 5110,
        "Fe2O3": 5250,
        "SiO2": 2320,
        "Al2O3": 3965,
        "CaO": 3340,
        "MgO": 3580,
        "C": 2267,
        "H2O": 1000,
        "Impurities": 3000,
    }

    # 5) Instantiate each step
    coarse_crusher = CoarseCrusher(xl_bin=">50 mm", large_bin="9.5-50 mm", medium_bin = "2.5-9.5 mm", fines_bin ="<2.5 mm", dust_loss_fraction=0.01)
    washer = Washer(removal_efficiency=1)

    dryer = Dryer(drying_efficiency=1.0)
    
    final_crusher = CrusherGrinder(xl_bin=">50 mm", large_bin="9.5-50 mm", medium_bin = "2.5-9.5 mm", fines_bin ="<2.5 mm",dust_loss_fraction=0.01)
    
    calcinator = Calcinator(
        heating_rate=10.0,
        final_temp=900.0,
        hold_time=60.0,
        gas_flow_rate=1.0,
        phase_densities=PHASE_DENSITIES
    )

    # 6) Build pipeline with recommended order
    pipeline = PretreatmentPipeline(
        throughput=throughput,
        coarse_crusher=coarse_crusher,
        washer=washer,
        final_crusher=final_crusher,
        dryer=dryer,
        calcinator=calcinator
    )

    # 7) Run pipeline
    final_ore = pipeline.run_pipeline(ore)

    # 8) Print final results
    print("\n=== Final Ore After Pretreatment ===")
    print(f"Total Mass: {final_ore.mass:.2f} kg")
    print(f"Temperature: {final_ore.temperature:.2f} °C")
    print(f"Density: {final_ore.density:.2f} kg/m^3")

    print("\nPhases (fraction of total mass):")
    for phase_name, fraction in final_ore.phases.items():
        print(f"  {phase_name}: {fraction*100:.2f}%")

    print("\nSize Distribution (fraction of total mass):")
    for bin_name, fraction in final_ore.size_distribution.items():
        print(f"  {bin_name}: {fraction:.2f}")

    # Print times
    total_time_hours = pipeline.total_time()
    print(f"\n--- Step Times ---")
    print(f"Coarse Crushing Time: {pipeline.time_coarse_crush:.2f} h")
    print(f"Washing Time:        {pipeline.time_washing:.2f} h")
    print(f"Final Crushing Time: {pipeline.time_final_crush:.2f} h")
    print(f"Drying Time:         {pipeline.time_drying:.2f} h")
    print(f"Calcination Time:    {pipeline.time_calcination:.2f} h")
    print(f"Total Time:          {total_time_hours:.2f} h")

    # If user selected mode '1', we can compute actual total run time
    # If user selected mode '2', we already used operation_time * throughput
    # so let's just display total_time anyway for reference.

    if calcinator:
        print(f"\nCO2 Released: {calcinator.total_co2_released:.2f} kg")

if __name__ == "__main__":
    main()
