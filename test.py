class S(dict):
    def __init__(self, *args, **kwargs):
        super().__init__()
        if kwargs:
            self.asdf = kwargs.get("asdf")

S(asdf=32, jk=23)

print(S)
