/**
 * Task Progress Card - Real-time task progress tracking
 */
import React, { useEffect, useState } from 'react'
import {
  Card,
  Progress,
  Space,
  Typography,
  Tag,
  Button,
  Tooltip,
  Alert,
  Statistic,
  Descriptions,
  Badge,
} from 'antd'
import {
  ClockCircleOutlined,
  CheckCircleOutlined,
  CloseCircleOutlined,
  SyncOutlined,
  StopOutlined,
  EyeOutlined,
  ThunderboltOutlined,
  FieldTimeOutlined,
} from '@ant-design/icons'
import { taskService } from '@/services/task.service'
import { toast } from '@/utils/notification'
import type { PipelineTask } from '@/types/task'
import dayjs from 'dayjs'
import duration from 'dayjs/plugin/duration'
import relativeTime from 'dayjs/plugin/relativeTime'

dayjs.extend(duration)
dayjs.extend(relativeTime)

const { Text, Title } = Typography

interface TaskProgressCardProps {
  taskId: string
  onComplete?: () => void
  onCancel?: () => void
  autoRefresh?: boolean
  refreshInterval?: number
}

export const TaskProgressCard: React.FC<TaskProgressCardProps> = ({
  taskId,
  onComplete,
  onCancel,
  autoRefresh = true,
  refreshInterval = 3000, // 3 seconds
}) => {
  const [task, setTask] = useState<PipelineTask | null>(null)
  const [loading, setLoading] = useState(true)
  const [elapsed, setElapsed] = useState(0)
  const [cancelling, setCancelling] = useState(false)

  // Load task details
  const loadTask = async () => {
    try {
      const data = await taskService.getTask(taskId)
      setTask(data)

      // Call onComplete callback if task completed
      if (
        data.status === 'completed' &&
        task?.status !== 'completed' &&
        onComplete
      ) {
        onComplete()
      }
    } catch (error: any) {
      console.error('Failed to load task:', error)
    } finally {
      setLoading(false)
    }
  }

  // Auto-refresh for running tasks
  useEffect(() => {
    loadTask()

    if (autoRefresh && task?.status === 'running') {
      const interval = setInterval(loadTask, refreshInterval)
      return () => clearInterval(interval)
    }
  }, [taskId, autoRefresh, task?.status])

  // Calculate elapsed time
  useEffect(() => {
    if (!task?.started_at) return

    const updateElapsed = () => {
      const start = dayjs(task.started_at)
      const now = dayjs()
      setElapsed(now.diff(start, 'second'))
    }

    updateElapsed()
    const interval = setInterval(updateElapsed, 1000)
    return () => clearInterval(interval)
  }, [task?.started_at])

  const handleCancel = async () => {
    if (!task) return

    setCancelling(true)
    try {
      await taskService.cancelTask(task.id)
      toast.success('Task cancelled successfully')
      await loadTask()
      if (onCancel) onCancel()
    } catch (error: any) {
      toast.error(`Failed to cancel task: ${error.message}`)
    } finally {
      setCancelling(false)
    }
  }

  // Status configuration
  const getStatusConfig = (status: string) => {
    const configs: Record<
      string,
      { color: string; icon: React.ReactNode; label: string; badge: string }
    > = {
      pending: {
        color: 'default',
        icon: <ClockCircleOutlined />,
        label: 'Pending',
        badge: 'default',
      },
      running: {
        color: 'processing',
        icon: <SyncOutlined spin />,
        label: 'Running',
        badge: 'processing',
      },
      completed: {
        color: 'success',
        icon: <CheckCircleOutlined />,
        label: 'Completed',
        badge: 'success',
      },
      failed: {
        color: 'error',
        icon: <CloseCircleOutlined />,
        label: 'Failed',
        badge: 'error',
      },
      cancelled: {
        color: 'default',
        icon: <StopOutlined />,
        label: 'Cancelled',
        badge: 'default',
      },
    }
    return (
      configs[status] || {
        color: 'default',
        icon: null,
        label: status,
        badge: 'default',
      }
    )
  }

  // Format elapsed time
  const formatElapsed = (seconds: number) => {
    const h = Math.floor(seconds / 3600)
    const m = Math.floor((seconds % 3600) / 60)
    const s = seconds % 60

    if (h > 0) return `${h}h ${m}m ${s}s`
    if (m > 0) return `${m}m ${s}s`
    return `${s}s`
  }

  if (loading || !task) {
    return <Card loading style={{ marginBottom: 16 }} />
  }

  const statusConfig = getStatusConfig(task.status)

  // Calculate estimated time remaining (mock calculation)
  const estimatedTotal = task.config?.estimated_time || 3600 // Default 1 hour
  const estimatedRemaining =
    task.status === 'running'
      ? Math.max(0, estimatedTotal - elapsed)
      : 0

  return (
    <Badge.Ribbon text={statusConfig.label} color={statusConfig.badge}>
      <Card
        size="small"
        title={
          <Space>
            {statusConfig.icon}
            <Text strong>{task.task_name}</Text>
          </Space>
        }
        extra={
          <Space>
            {task.status === 'running' && (
              <Button
                size="small"
                danger
                icon={<StopOutlined />}
                onClick={handleCancel}
                loading={cancelling}
              >
                Cancel
              </Button>
            )}
            <Button size="small" icon={<EyeOutlined />}>
              View Details
            </Button>
          </Space>
        }
        style={{ marginBottom: 16 }}
      >
        <Space direction="vertical" size="middle" style={{ width: '100%' }}>
          {/* Progress Bar */}
          {task.status === 'running' && (
            <div>
              <div style={{ marginBottom: 8, display: 'flex', justifyContent: 'space-between' }}>
                <Text type="secondary">Progress</Text>
                <Text type="secondary">{Math.round(task.progress || 0)}%</Text>
              </div>
              <Progress
                percent={Math.round(task.progress || 0)}
                status={task.status === 'failed' ? 'exception' : 'active'}
                strokeColor={{
                  '0%': '#108ee9',
                  '100%': '#87d068',
                }}
              />
            </div>
          )}

          {/* Time Statistics */}
          {task.status === 'running' && (
            <Space size="large" wrap>
              <Statistic
                title="Elapsed"
                value={formatElapsed(elapsed)}
                prefix={<FieldTimeOutlined />}
                valueStyle={{ fontSize: 14 }}
              />
              {estimatedRemaining > 0 && (
                <Statistic
                  title="Estimated Remaining"
                  value={formatElapsed(estimatedRemaining)}
                  prefix={<ClockCircleOutlined />}
                  valueStyle={{ fontSize: 14 }}
                />
              )}
            </Space>
          )}

          {/* Completed */}
          {task.status === 'completed' && task.completed_at && (
            <Alert
              message="Task Completed Successfully"
              description={`Finished at ${dayjs(task.completed_at).format(
                'YYYY-MM-DD HH:mm:ss'
              )}`}
              type="success"
              showIcon
              icon={<CheckCircleOutlined />}
            />
          )}

          {/* Failed */}
          {task.status === 'failed' && task.error_message && (
            <Alert
              message="Task Failed"
              description={task.error_message}
              type="error"
              showIcon
              icon={<CloseCircleOutlined />}
            />
          )}

          {/* Task Details */}
          <Descriptions size="small" column={2} bordered>
            <Descriptions.Item label="Task Type">
              <Tag icon={<ThunderboltOutlined />} color="blue">
                {task.task_type || 'N/A'}
              </Tag>
            </Descriptions.Item>
            <Descriptions.Item label="Started">
              {task.started_at
                ? dayjs(task.started_at).format('YYYY-MM-DD HH:mm:ss')
                : 'Not started'}
            </Descriptions.Item>
            {task.celery_task_id && (
              <Descriptions.Item label="Celery ID" span={2}>
                <Tooltip title={task.celery_task_id}>
                  <Text
                    copyable={{ text: task.celery_task_id }}
                    style={{ fontSize: 12, fontFamily: 'monospace' }}
                  >
                    {task.celery_task_id.slice(0, 16)}...
                  </Text>
                </Tooltip>
              </Descriptions.Item>
            )}
          </Descriptions>

          {/* Live Updates Indicator */}
          {task.status === 'running' && autoRefresh && (
            <div style={{ textAlign: 'center' }}>
              <Text type="secondary" style={{ fontSize: 12 }}>
                <SyncOutlined spin style={{ marginRight: 4 }} />
                Live updates every {refreshInterval / 1000}s
              </Text>
            </div>
          )}
        </Space>
      </Card>
    </Badge.Ribbon>
  )
}

export default TaskProgressCard
