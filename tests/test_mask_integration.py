from lagkinematic.integration import EulerIntegrator, ParticleState, ParticlePairState
from lagkinematic.geometry import LonLatDomain

def test_integration_mask():
    domain = LonLatDomain(lon_min=0, lon_max=360, lat_min=-90, lat_max=90)

    class FakeMask:
        def sample_mask(self, lon, lat, depth):
            return 0.0 if lat > 37 else 1.0

    class FakeSampler:
        def sample(self, lon, lat, depth, t):
            return 0.1, 0.0, 0.0  # 0.1 m/s verso Est

    mask = FakeMask()
    sampler = FakeSampler()

    integrator = EulerIntegrator(
        sampler=sampler,
        domain=domain,
        dt=3600.0,
        mask=mask,
        mask_threshold=0.5,
        subgrid=None,
    )

    p1 = ParticleState(id=1, lon=12.0, lat=36.0, depth=0.0)  # mare
    p2 = ParticleState(id=2, lon=12.0, lat=38.0, depth=0.0)  # terra

    pair = ParticlePairState(id=0, p1=p1, p2=p2, t=0.0)
    pairs = [pair]

    integrator.apply_initial_mask(pairs)

    print("p1 alive:", p1.alive, "p2 alive:", p2.alive)

    assert p1.alive is False
    assert p2.alive is False

    integrator.step_pair(pair, 0.0)

    print("p1 lon after:", p1.lon, "alive:", p1.alive)
    print("p2 lon after:", p2.lon, "alive:", p2.alive)

    assert p1.alive is False
    assert p2.alive is False


def test_integration_mask_off_mode_skips_beaching():
    domain = LonLatDomain(lon_min=0, lon_max=360, lat_min=-90, lat_max=90)

    class FakeMask:
        def sample_mask(self, lon, lat, depth):
            return 0.0

    class FakeSampler:
        def sample(self, lon, lat, depth, t):
            return 0.0, 0.1, 0.0  # 0.1 m/s verso Nord

    mask = FakeMask()
    sampler = FakeSampler()

    integrator = EulerIntegrator(
        sampler=sampler,
        domain=domain,
        dt=3600.0,
        mask=mask,
        mask_threshold=0.5,
        subgrid=None,
        beaching_mode="off",
    )

    p1 = ParticleState(id=1, lon=12.0, lat=36.0, depth=0.0)
    p2 = ParticleState(id=2, lon=12.0, lat=36.0, depth=0.0)

    pair = ParticlePairState(id=0, p1=p1, p2=p2, t=0.0)
    pairs = [pair]

    integrator.apply_initial_mask(pairs)

    assert p1.alive is True
    assert p2.alive is True

    integrator.step_pair(pair, 0.0)

    assert p1.alive is True
    assert p2.alive is True
    assert p1.lat > 36.0

if __name__ == "__main__":
    test_integration_mask()
