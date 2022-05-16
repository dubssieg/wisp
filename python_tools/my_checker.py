# effectue des tests de types d'entrée/sortie de fonctions et de classes

from inspect import signature, _empty


class TCError(BaseException):
    """
    Erreur custom permettant d'afficher un message d'erreur.
    """

    def __init__(self, msg: str) -> None:
        self.__message = msg
        super().__init__(self.__message)


def my_types_checker(func):
    """
    Attendu pour être utilisé comme décorateur.
    ---
    Effectue une vérifiaction sur les types explicités dans une signature de fonction ou de classe décorée.
    Si le type ne correspond pas au type signifié, on lève une erreur
    Si le type n'est pas spécifié, on ignore l'argument.
    Fonctionne à la fois avec les args et les kwargs.
    Permet également de vérifier le type de retour, ignoré également si il n'est pas indiqué.
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
        for input_type, annotation in comp:
            if annotation is not _empty and input_type != annotation:
                raise TCError(
                    f"Erreur de la fonction {func.__name__} : le type entré {input_type} ne correspond pas au type attendu {annotation}")

        # execution de la fonction décorée
        retour = func(*args, **kwargs)

        # comparaison du type de retour
        if ret_annotation is not _empty and ret_annotation != None and type(retour) != ret_annotation:
            raise TCError(
                f"Erreur de la fonction {func.__name__} : le type de retour {type(retour)} ne correspond pas au type attendu {ret_annotation}")

        return retour
    return wrapper


def my_class_checker(cls):
    """
    Application du décorateur à l'entèreté d'une classe.
    A ne pas utiliser sur une fonction seul. 
    """
    for attr in cls.__dict__:
        if callable(getattr(cls, attr)):
            setattr(cls, attr, my_types_checker(getattr(cls, attr)))
    return cls
