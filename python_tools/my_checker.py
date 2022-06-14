# IO tests on classes level

from inspect import signature, _empty
from typing import Callable


def my_types_checker(func: Callable):
    """Checks the types of annotations of a func and raises error if not corresponding.

    Args:
        func (Callable): a func
    """
    def wrapper(*args, **kwargs):
        """
        Corps de test
        """

        # formation des listes de comparaison
        arguments = args + tuple(kwargs.keys())
        args_types = [type(arg) for arg in arguments]
        signature_func = signature(func).parameters
        ret_annotation = signature(func).return_annotation
        annotations = [
            signature_func[elt].annotation for elt in signature_func]

        # comparaison des annotations de signature et des types des arguments
        comp = zip(args_types, annotations)
        for i, (input_type, annotation) in enumerate(comp):
            if annotation is not _empty and input_type != annotation:
                try:
                    annotation(arguments[i])
                except Exception as exc:
                    raise BaseException(
                        f"Error in {func.__name__} : input type {input_type} incompatible with {annotation}") from exc

        # execution de la fonction décorée
        retour = func(*args, **kwargs)

        # comparaison du type de retour
        if ret_annotation is not _empty and ret_annotation is not None and isinstance(retour, ret_annotation):
            raise BaseException(
                f"Error in {func.__name__} : return type {type(retour)} does not match {ret_annotation}")

        return retour
    return wrapper


def my_class_checker(arg: Callable):
    def wrapper(cls):
        """Application of a decorator on whole class

        arg is a callable, decorator to apply

        Returns:
            cls : decorated classes
        """
        for attr in cls.__dict__:
            if callable(getattr(cls, attr)):
                setattr(cls, attr, arg(getattr(cls, attr)))
        return cls
    return wrapper
