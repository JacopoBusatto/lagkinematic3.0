from lagkinematic.integration import EulerIntegrator, ParticleState, ParticlePairState
from lagkinematic.geometry import LonLatDomain


class FakeBottom:
    def __init__(self, depth: float = 5.0):
        self.depth = depth

    def sample_bottom(self, lon, lat):
        return float(self.depth)


class FakeSampler:
    def __init__(self, w: float):
        self.w = w

    def sample(self, lon, lat, depth, t):
        return 0.0, 0.0, self.w


class FakeMask:
    def __init__(self, value: float):
        self.value = value

    def sample_mask(self, lon, lat, depth):
        return self.value


def test_bottom_kill_mode():
    domain = LonLatDomain(lon_min=0, lon_max=10, lat_min=0, lat_max=10)
    sampler = FakeSampler(w=2.0)
    bottom = FakeBottom(depth=5.0)

    p1 = ParticleState(id=1, lon=1.0, lat=1.0, depth=4.0)
    p2 = ParticleState(id=2, lon=1.0, lat=1.0, depth=1.0)
    pair = ParticlePairState(id=0, p1=p1, p2=p2)

    integrator = EulerIntegrator(
        sampler=sampler,
        domain=domain,
        dt=1.0,
        bottom=bottom,
        bottom_mode="kill",
    )

    integrator.step_pair(pair, t=0.0)

    assert p1.alive is False
    assert p2.alive is False


def test_bottom_bounce_mode():
    domain = LonLatDomain(lon_min=0, lon_max=10, lat_min=0, lat_max=10)
    sampler = FakeSampler(w=2.0)
    bottom = FakeBottom(depth=5.0)

    p1 = ParticleState(id=1, lon=1.0, lat=1.0, depth=4.0)
    p2 = ParticleState(id=2, lon=1.0, lat=1.0, depth=1.0)
    pair = ParticlePairState(id=0, p1=p1, p2=p2)

    integrator = EulerIntegrator(
        sampler=sampler,
        domain=domain,
        dt=1.0,
        bottom=bottom,
        bottom_mode="bounce",
    )

    integrator.step_pair(pair, t=0.0)

    assert p1.alive is True
    assert p2.alive is True
    assert p1.depth <= bottom.depth
    assert p1.depth < 5.1


def test_bottom_off_mode_allows_exceeding_depth():
    domain = LonLatDomain(lon_min=0, lon_max=10, lat_min=0, lat_max=10)
    sampler = FakeSampler(w=3.0)
    bottom = FakeBottom(depth=2.0)

    p1 = ParticleState(id=1, lon=1.0, lat=1.0, depth=1.5)
    p2 = ParticleState(id=2, lon=1.0, lat=1.0, depth=1.5)
    pair = ParticlePairState(id=0, p1=p1, p2=p2)

    integrator = EulerIntegrator(
        sampler=sampler,
        domain=domain,
        dt=1.0,
        bottom=bottom,
        bottom_mode="off",
    )

    integrator.step_pair(pair, t=0.0)

    assert p1.alive is True
    assert p2.alive is True
    assert p1.depth > bottom.depth


def test_apply_initial_mask_kills_on_land_even_with_bounce():
    domain = LonLatDomain(lon_min=0, lon_max=10, lat_min=0, lat_max=10)
    sampler = FakeSampler(w=0.0)
    mask = FakeMask(value=0.0)  # terra

    p1 = ParticleState(id=1, lon=1.0, lat=1.0, depth=0.0)
    p2 = ParticleState(id=2, lon=1.0, lat=1.0, depth=0.0)
    pair = ParticlePairState(id=0, p1=p1, p2=p2)
    pairs = [pair]

    integrator = EulerIntegrator(
        sampler=sampler,
        domain=domain,
        dt=1.0,
        mask=mask,
        mask_threshold=0.5,
        beaching_mode="bounce",
        bottom_mode="bounce",
    )

    integrator.apply_initial_mask(pairs)

    assert p1.alive is False
    assert p2.alive is False
