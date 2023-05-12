from main import sum
import pytest


def test_sum():
    assert sum(1, 2) == 3


def test_sum_diminuir():
    assert sum(1, 2)-2 == 1
