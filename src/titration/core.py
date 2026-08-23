from .models import TitrationInput, TitrationResult
from .solver import run_solver


def run_titration(inp: TitrationInput) -> TitrationResult:
    return run_solver(inp)
