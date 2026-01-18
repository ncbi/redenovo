import sys
from timeit import default_timer as timer
import numpy as np
import tensorflow as tf

# get logger instance
import logging
_logger = logging.getLogger(__name__)


class Optimizer(object):

    # flattens heterogeneous list
    def flat2gen(self, alist):
        for item in alist:
            if isinstance(item, list):
                for subitem in item: yield subitem
            else:
                yield item

    def __init__(self, data):
        _logger.debug("Creating instance of Optimizer")

        self.data = data
        self.args = data.args

        # Get a tf.Variable representation of all data required for the optimization
        self.M = tf.Variable(initial_value=tf.convert_to_tensor(value=data.get_m()),
                            trainable=False,
                            name='M')

        self.P_fixed, self.P_inferred = self.make_p()
        _logger.debug("Finished make p, trying to make A")
        self.A = tf.Variable(initial_value=tf.convert_to_tensor(value=data.get_a()),
                            trainable=True,
                            constraint=lambda x: tf.clip_by_value(x,0,np.infty),
                            name='A')

        _logger.debug("Completed tf.Variable initialization")

        # Define which tensors should be optimized
        self.trainable_variables = list(filter(lambda x: x is not None and x.trainable, self.flat2gen([self.P_inferred, self.A])))

        _logger.debug(f"Trainable variables are {[x.name for x in self.trainable_variables]}")
        
        # Global spaceholder for the chosen optimizer
        self.optimizer = None

        # store performance values of interest here
        self.logs = {
                    'losses'      : tf.Variable([], dtype=tf.float64, shape=tf.TensorShape(None)),
                    'best_losses' : tf.Variable([], dtype=tf.float64, shape=tf.TensorShape(None)),
                    'epochs'      : tf.Variable([], dtype=tf.int32,   shape=tf.TensorShape(None)),
                    'steps'       : tf.Variable([], dtype=tf.float32, shape=tf.TensorShape(None))
                    }

        _logger.debug("Successfully initialized optimizer.")

    def make_p(self):
        """ Creates a Variable of P with corresponding contraints

            Array corrsponding to the fixed and inferred portion of P
        """

        _logger.debug("make_p")

        # get the fixed and inferred portion and combine
        p = self.data.get_p()

        p_fixed = p[0]  # this is either None or not
        if p_fixed is not None:  # Convert to Variable
            p_fixed = tf.Variable(
                                    initial_value=tf.convert_to_tensor(value=p[0]),
                                    trainable=False,
                                    name='P_fixed'
                                )


        p_inferred = p[1]
        if p_inferred is not None:
            # constraint: probabilities must sum to one along the mutational signature axis (1)
            clip_min = 1e-10
            clip_max = np.infty
            sum_one_axis = 1
            constraint = lambda x: tf.clip_by_value(x, clip_min, clip_max) / tf.reduce_sum(tf.clip_by_value(x, clip_min,clip_max), axis=sum_one_axis, keepdims=True)
            p_inferred = tf.Variable(
                                initial_value=tf.convert_to_tensor(value=p[1]),
                                trainable=True,
                                constraint=constraint,
                                name='P_inferred')

        return [p_fixed, p_inferred]

    @tf.function
    def compute_loss(self):
        """ Returns an objective funtion corresponding to the chosen settings in args
        """
        P = tf.concat(list(filter(lambda x: x is not None, [self.P_fixed, self.P_inferred])), axis=0)
        Mhat = tf.matmul(self.A, P)
        loss_value = tf.norm(
                    tensor=self.M-Mhat,
                    ord='euclidean',
                    name='frobenius_norm'
                    )
        return loss_value

    
    @tf.function
    def train_step(self):
        """ Wrap loss and gradient computation in a tf.function to speed up computation.
            This should lead to at least 5 times faster running times as things are now not
            executed eagerly but compiled into the compute graph.
        """

        # Compute loss using a GradientTape
        with tf.GradientTape() as tape:
            loss = self.compute_loss()

        # Compute gradients
        gradients = tape.gradient(loss, self.trainable_variables)

        # Apply gradients
        self.optimizer.apply_gradients(zip(gradients, self.trainable_variables))

        return loss
        
    def optimize(self):
        """
        Contains logic to perform the full optimization
        """

        _logger.debug("Model definition complete, commencing optimization")
        _logger.debug("Optimization Strategy:")
        for epoch, stepsize, optimizer in zip(self.args.epochs, self.args.stepsizes, self.args.optimizers):
            _logger.debug(f"{epoch} iterations with step size {stepsize}; optimizer {optimizer}")

        # Perform optimization
        t0=timer()
        for epoch, stepsize, optimizer in zip(self.args.epochs, self.args.stepsizes, self.args.optimizers):
            _logger.debug(f"### New optimization: iterations {epoch} step size {stepsize} optimizer {optimizer}")
            #tf.print(f"### New optimization: iterations {epoch} step size {stepsize} optimizer {optimizer}")
            self.optimizeModel(self.compute_loss, epoch, stepsize, optimizer)
            _logger.debug(f"Best loss {self.compute_loss()}")
        t1 = timer()

        # Wall seconds elapsed (floating point)
        _logger.debug(f"Optimization complete. Time elapsed: {t1 - t0} seconds")
        return self.compute_loss()
        

    def optimizeModel(self, loss, iters, stepsize, type):
        """
        Performs one optimization bracket according to the provided options

        Args:
            loss: identifier of which loss function to use
            iters: number of iterations (epochs) to perform
            stepsize: initial stepsize for the optimizer
            type: the optimizer type to use
        """

        _logger.debug("Preparing optimization")

        # Convert everything into tensorfow equivalents so we can run it on the GPU

        # print directives
        update_steps = tf.constant(self.args.optimizer_user_update_steps)
        log_update_steps = tf.constant(self.args.optimizer_log_update_steps)

        # optimization details
        epoch = tf.Variable(0)
        num_it = tf.constant(iters)
        cont = tf.Variable(True)

        # current best loss
        best_loss = tf.Variable(loss()) #best_loss = tf.Variable(tf.cast(loss(), tf.float32))

        # best parameters so far. Since these change depending on the user input,
        # use an array to store them. They should always correspond to the
        # content and order of self.trainable_variables
        parameters_best = [tf.Variable(tf.identity(x), shape=x.shape) for x in self.trainable_variables]

        # Define a training operation for tensforflow, this can be exchanged with other optimizers if desired
        # adadelta,adagrad,adam,adamax,nadam,rmaprop,sgd
        stepsize = float(stepsize)
        if type == "adadelta":
            self.optimizer = tf.optimizers.Adadelta(stepsize)
        elif type == "adagrad":
            self.optimizer = tf.optimizers.Adagrad(stepsize)
        elif type == "adam":
            self.optimizer = tf.optimizers.Adam(stepsize)
        elif type == "adamax":
            self.optimizer = tf.optimizers.Adamax(stepsize)
        elif type == "nadam":
            self.optimizer = tf.optimizers.Nadam(stepsize)
        elif type == "rmaprop":
            self.optimizer = tf.optimizers.RMSprop(stepsize)
        elif type == "sgd":
            self.optimizer = tf.optimizers.SDG(stepsize)
        elif type == _:
            self.optimizer = tf.optimizers.Adam(stepsize)

        _logger.debug(f"Optimizer is {self.optimizer}")

        # see https://www.tensorflow.org/guide/keras/writing_a_training_loop_from_scratch for more options on custom training loops
        while epoch < num_it and cont:

            loss = self.train_step()

            # store the best weights so far, these are the ones we will return
            if loss < best_loss:
                best_loss.assign(loss)
                for i, x in enumerate(self.trainable_variables):
                    parameters_best[i].assign(x)

            # update logs?
            if epoch % log_update_steps == 0:
                self.logs['losses'] = tf.concat([self.logs['losses'], tf.Variable([loss])], axis=0)
                self.logs['best_losses'] = tf.concat([self.logs['best_losses'], tf.Variable([best_loss])], axis=0)
                self.logs['epochs'] = tf.concat([self.logs['epochs'], tf.Variable([epoch])], axis=0)
                self.logs['steps'] = tf.concat([self.logs['steps'], tf.Variable([stepsize])], axis=0)

            # print update?
            #tf.cond(tf.equal(epoch % update_steps, 0),
            #        lambda: [   #True
            #                    tf.print("Epoch", epoch, "loss", loss, "best loss", best_loss, output_stream=sys.stdout)
            #                ],
            #        lambda: [   #False
            #                    tf.no_op()  # just a placeholder
            #                ])

            epoch.assign(epoch + 1)

        # Save best results at the end of the current iteration
        for x, y in zip(parameters_best, self.trainable_variables):
            y.assign(x)

        #tf.print("Optimization bracket completed")
        _logger.debug("Optimization bracket completed")
        
        
    def store(self):
        """
        Stores the optimized tensors back in the resources class. Takes care of
        removing the extra dimensions

        """
        _logger.debug(f" Storing optimized parameters")

        self.data.set_p([
                        None if self.P_fixed is None else self.P_fixed.numpy(),
                        None if self.P_inferred is None else self.P_inferred.numpy()
                        ])

        self.data.set_a(self.A.numpy())

        self.data.set_logs({log: self.logs[log].numpy() for log in self.logs})
