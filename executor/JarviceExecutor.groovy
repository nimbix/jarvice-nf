package nextflow.executor
import java.nio.file.Path

import groovy.util.logging.Slf4j
import nextflow.processor.TaskRun

@Slf4j
class JarviceExecutor extends AbstractGridExecutor {

    // Jarvice specific cluster options
    private String jsonClusterOpts;
    private String sPipelineID;

    @Override
    protected String getHeaderToken() { return '#JARVICE' }

    @Override
    protected List<String> getDirectives(TaskRun task, List<String> result)
    {
        /*
        // Handle common executor directives (long options used for clarity)
        result << '--work-dir' << quote(task.workDir)
        result << '--job-name' << getJobNameFor(task)
        result << '--cmd-log' << quote(task.workDir.resolve(TaskRun.CMD_LOG))

        // Jarvice specific directives
        if( task.config.container ) {
            result << '--app ' + task.config.container.toString()
        }

        if( task.config.machine ) {
            result << '--machine ' + task.config.machine.toString()
        }

        if( task.config.label ) {
            result << '--job-label ' + task.config.label.toString()
        }

        if( task.config.time ) {
            result << '--wall-time ' + task.config.time.toString()
        }

        // Jarvice related stuff
        if( task.config.user ) {
            result << '--user ' + task.config.user.toString()
        }

        if( task.config.apikey ) {
            result << '--apikey ' + task.config.apikey.toString()
        }

        if( task.config.apiserver ) {
            result << '--apiserver ' + task.config.apiserver.toString()
        }

        log.debug "getDirectives $result"
        */

        // -- at the end append the command script wrapped file name
        /*if( task.config.clusterOptions ) {
            result << task.config.clusterOptions.toString() << ''
        }*/

        return result
    }

    @Override
    List<String> getSubmitCommandLine(TaskRun task, Path scriptFile)
    {
        log.debug "task $task"
        sPipelineID = task.processor.getSession().runName;

        List<String> result = ['jarvice-nxf-submit']
        if(task.config.machineType) result << '--machine' << task.config.machineType.toString();
        if(task.config.cpus) result << '--nodes' << task.config.cpus.toString();
        result << '--label' << sPipelineID;

        // Save cluster options so the queue status command can authenticate
        if(task.config.clusterOptions)
        {
            jsonClusterOpts = task.config.clusterOptions.toString()
            result << '--opts' << jsonClusterOpts
        };

        result << scriptFile.getName()
        return result
    }

    @Override
    def parseJobId(String text)
    {
        log.debug "parseJobId $text"
        // TODO parse submission errors and throw exception if submit failed
        return text.trim()
    }

    @Override
    protected List<String> getKillCommand() {
        log.debug "getKillCommand"
        ['jarvice-nxf-kill']
    }

    @Override
    protected List<String> queueStatusCommand(Object queue)
    {
        log.debug "queueStatusCommand $queue"
        ['jarvice-nxf-qstatus', sPipelineID]
    }

    /*
    Job status and substatus constants reference
        SUBSTATUS_NONE: 0,
        SUBSTATUS_LIMITS: 1,
        SUBSTATUS_TERMINATED_OUTSIDE: 2,
        SUBSTATUS_LICENSE_FEATURES: 3,
        SUBSTATUS_SUSPENDED_USER: 4,
        SUBSTATUS_SUSPENDED_TEMP_LM: 5,

        SUBMITTED             - Job submitted successfully
        PROCESSING STARTING   - Job is running
        COMPLETED             - Job ended with exit status 0
        COMPLETED WITH ERROR  - Job ended with non-zero exit status
        TERMINATED            - Job ended with SIGTERM (user shut down the job)
        CANCELED              - Job ended with other other signal
        EXEMPT                - Job will be exempt from user billing
        SEQUENTIALLY QUEUED   - Job will be queued after other jobs (future use)
     */

    // Map of jarvice status:substatus to Nextflow
    // Note that jarvice-nxf-qstatus replaces spaces with _ for the job status strings
    static private Map STATUS_MAP = [
            'SUBMITTED:0': QueueStatus.PENDING,     // Waiting for job to start
            'SUBMITTED:1': QueueStatus.HOLD,        // Held up waiting for user limits
            'SUBMITTED:3': QueueStatus.HOLD,        // Held up waiting for license features

            'PROCESSING_STARTING:0': QueueStatus.RUNNING,   // Running
            'PROCESSING_STARTING:4': QueueStatus.HOLD,      // Suspended by user
            'PROCESSING_STARTING:5': QueueStatus.HOLD,      // Suspended by License manager

            'COMPLETED:0': QueueStatus.DONE,
            'COMPLETED_WITH_ERROR:0': QueueStatus.ERROR,

            'CANCELED:0': QueueStatus.ERROR,        // Killed by Jarvice before job started
            'CANCELED:2': QueueStatus.ERROR,        // Killed before job started outside of Jarvice

            'TERMINATED:0': QueueStatus.ERROR,      // Killed by Jarvice after job started
            'TERMINATED:2': QueueStatus.ERROR,      // Killed outside of Jarvice after job started

            'EXEMPT:0': QueueStatus.DONE,                       // Completed job marked as exempt from payment manually
            'SEQUENTIALLY_QUEUED:0': QueueStatus.PENDING,       // Not used currently
    ]

    @Override
    protected Map<String, QueueStatus> parseQueueStatus(String text)
    {
        def result = [:]
        text.eachLine {
            String line -> def cols = line.split(/\s+/)
            if( cols.size() == 4 )
            {
                result.put( cols[0], STATUS_MAP.get(cols[1]) )
            }
            else
            {
                log.debug "[JARVICE] invalid status line: `$line`"
            }
        }

        return result
    }
}
