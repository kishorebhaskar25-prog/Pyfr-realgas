def call_(obj, name_, **kwargs):
    keys = list(kwargs)
    keys[0] = keys[0][0].upper() + keys[0][1:]
    name_ = str(name_)
    meth = f"{name_}{'_'.join(str(k) for k in keys)}_"

    return getattr(obj, meth)(*kwargs.values())


def init_(cls, **kwargs):
    return call_(cls.alloc(), 'initWith', **kwargs)
