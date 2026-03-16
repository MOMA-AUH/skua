from skua import hello


def test_hello_default() -> None:
    assert hello() == "Hello, world!"


def test_hello_custom_name() -> None:
    assert hello("Skua") == "Hello, Skua!"
