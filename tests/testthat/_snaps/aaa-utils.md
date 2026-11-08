# enzyme list processing rejects mixed concrete and abstract enzymes

    Code
      .process_enzymes_arg(c("MAN2A1", "ManII"))
    Condition
      Error in `.enzymes_from_list()`:
      ! Concrete and abstract enzymes cannot be mixed in `enzymes`.

---

    Code
      .process_enzymes_arg(list(enzyme("MAN2A1"), enzyme("ManII")))
    Condition
      Error in `.enzymes_from_list()`:
      ! Concrete and abstract enzymes cannot be mixed in `enzymes`.
