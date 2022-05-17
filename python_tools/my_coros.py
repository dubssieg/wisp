from queue import Queue
from threading import Thread, Semaphore, Lock, currentThread
from concurrent.futures import ThreadPoolExecutor
from typing import Callable


def my_coro(func, args: list):
    """
    Shueldues a list of coroutines
    * func is a single func which will be executed len(args) times
    * args is a list of serial args for func (eg list of list)
    """
    handlers = [Queue() for _ in range(len(args))]
    funcs = [my_thread(handlers[i])(func)(*args[i]) for i in range(len(args))]
    return [hdl.get() for hdl in handlers]


def my_thread(resultQueue=None):
    def wrapper(function: Callable):
        def pw(*args, **kwargs):
            def process(*args, **kwargs):
                ret = function(*args, **kwargs)
                if resultQueue:
                    resultQueue.put(ret)
            thread = Thread(target=process, args=args, kwargs=kwargs)
            thread.setDaemon(True)
            thread.start()
            return process
        return pw
    return wrapper


def parallelize(func: Callable, list_args: list[tuple]) -> list:
    list_returns: list = []
    pool = ThreadPool()
    sem = Semaphore(8)
    for i in range(len(list_args)):
        t = Thread(target=launcher, name=f"thread_{i}", args=(
            sem, pool, func, list_args[i]))
        t.start()
    return list_returns


class ThreadPool(object):
    def __init__(self):
        super(ThreadPool, self).__init__()
        self.active = []
        self.lock = Lock()

    def makeActive(self, name):
        with self.lock:
            self.active.append(name)

    def makeInactive(self, name):
        with self.lock:
            self.active.remove(name)


def launcher(semaphore: Semaphore, pool_exec: ThreadPool, func: Callable, args):
    with semaphore:
        name = currentThread().getName()
        pool_exec.makeActive(name)
        ret = func(*args)
        pool_exec.makeInactive(name)
        return ret


def my_futures_collector(func: Callable, argslist: list, num_processes: int) -> list:
    """
    Spawns len(arglist) instances of func and executes them at num_processes instances at time.

    * func : a function
    * argslist (list): a list of tuples, arguments of each func
    * num_processes (int) : max number of concurrent instances
    """
    with ThreadPoolExecutor(max_workers=num_processes) as executor:
        futures = [executor.submit(func, *args) for args in argslist]
    return [f.result() for f in futures]
