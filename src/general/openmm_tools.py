from openmm.openmm import System, VerletIntegrator, Context
import openmm.unit as units


def create_dummy_context(platform, properties=None):
    system = System()
    system.addParticle(1.0)
    step_size_Ps = 0.001 * units.picosecond
    integrator = VerletIntegrator(step_size_Ps)
    if properties is None:
        context = Context(system, integrator, platform)
    else:
        context = Context(system, integrator, platform, properties)
    return platform, context