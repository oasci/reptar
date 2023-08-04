import time
from reptar import File
from reptar.calculators import Data, Driver
from reptar.calculators.psi4_workers import psi4_worker

rfile = File("1h2o.zarr", mode="a")

tasks = ["E", "G"]

data = Data(rfile)
data.prepare_tasks(
    tasks=tasks,
    source_key="/",
    source_labels={"Z": "atomic_numbers", "R": "geometry"},
    dest_key="/",
    dest_labels={"E": "energy_ele_df.mp2.def2tzvppd", "G": "grads_df.mp2.def2tzvppd"},
)
data.validate(tasks=tasks)

driver = Driver(
    use_ray=False, n_workers=1, n_cpus_per_worker=2, chunk_size=1, ray_address="auto"
)

worker_kwargs = {
    "charge": 0,
    "mult": 1,
    "method": "mp2",
    "threads": 2,
    "mem": "4 GB",
    "options": {
        "reference": "rhf",
        "scf_type": "df",
        "mp2_type": "df",
        "e_convergence": 10,
        "d_convergence": 10,
        "basis": "def2-tzvppd",
        "df_basis_scf": "def2-universal-jkfit",
        "df_basis_mp2": "def2-tzvppd-ri",
        "freeze_core": True,
    },
}

t_start = time.time()
data = driver.run(
    psi4_worker,
    worker_kwargs,
    data,
    tasks,
    start_slice=None,
    end_slice=5,
)
t_stop = time.time()
print("\n\n")
print(f"Energies: {data.E[:5]}")
print("Gradients:")
print(data.G[:5])
print(f"Calculations took {t_stop-t_start:.1f} seconds")
