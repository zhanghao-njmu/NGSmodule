/**
 * Task List Page - Task monitoring with real-time updates
 * Modernized with animations and enhanced UI components
 */
import type React from 'react'
import { useEffect, useState } from 'react'
import { useNavigate } from 'react-router-dom'
import { Button, Tag, Progress, Space, Select, Tooltip, Typography, Drawer, Alert, Spin } from 'antd'
import {
  PlusOutlined,
  SyncOutlined,
  CheckCircleOutlined,
  CloseCircleOutlined,
  StopOutlined,
  EyeOutlined,
  ThunderboltOutlined,
  FileTextOutlined,
  CopyOutlined,
} from '@ant-design/icons'
import type { ColumnsType } from 'antd/es/table'
import { useTaskStore } from '../../store/taskStore'
import { useProjectStore } from '../../store/projectStore'
import { websocketService } from '../../services/websocket.service'
import { taskService } from '../../services/task.service'
import {
  PageHeader,
  DataTable,
  StatisticCard,
  StatusTag,
  PageSkeleton,
  FadeIn,
  EnhancedEmptyState,
  StaggeredList,
} from '../../components/common'
import { confirm } from '../../components/common/ConfirmDialog'
import { toast } from '../../utils/notification'
import type { StatisticItem } from '../../components/common'
import type { Task, TaskStatus } from '../../types/task'
import dayjs from 'dayjs'

const { Option } = Select
const { Title, Text } = Typography

export const TaskList: React.FC = () => {
  const navigate = useNavigate()
  const { tasks, stats, loading, fetchTasks, fetchStats, cancelTask, handleWebSocketMessage } = useTaskStore()
  const { items, fetchItems } = useProjectStore()
  const [selectedProject, setSelectedProject] = useState<string>('')
  const [initialLoad, setInitialLoad] = useState(true)

  // Log viewer state
  const [logDrawerVisible, setLogDrawerVisible] = useState(false)
  const [currentTask, setCurrentTask] = useState<Task | null>(null)
  const [taskLogs, setTaskLogs] = useState<string>('')
  const [loadingLogs, setLoadingLogs] = useState(false)

  useEffect(() => {
    const loadData = async () => {
      await fetchItems()
      await fetchStats()
      await fetchTasks()
      setInitialLoad(false)
    }

    loadData()

    // Setup WebSocket for real-time updates
    const token = localStorage.getItem('auth_token')
    if (token) {
      websocketService.connect(token)
      websocketService.addMessageHandler(handleWebSocketMessage)
    }

    return () => {
      websocketService.disconnect()
    }
  }, [])

  useEffect(() => {
    if (!initialLoad) {
      if (selectedProject) {
        fetchTasks({ project_id: selectedProject })
      } else {
        fetchTasks()
      }
    }
  }, [selectedProject, initialLoad])

  const handleCancelTask = (task: Task) => {
    confirm({
      title: '取消任务',
      content: `您确定要取消任务 "${task.task_name}" 吗？`,
      okText: '取消任务',
      cancelText: '保持运行',
      okButtonProps: { danger: true },
      onConfirm: async () => {
        const loadingToast = toast.loading('取消中...')
        try {
          await cancelTask(task.id)
          loadingToast()
          toast.success('任务已取消')
          fetchTasks()
          fetchStats()
        } catch (error) {
          loadingToast()
          toast.error('取消失败，请重试')
        }
      },
    })
  }

  const handleViewLogs = async (task: Task) => {
    setCurrentTask(task)
    setLogDrawerVisible(true)
    setTaskLogs('')
    setLoadingLogs(true)

    try {
      const response = await taskService.getTaskLogs(task.id)
      setTaskLogs(response.log_content || 'No logs available')
    } catch (error: any) {
      setTaskLogs(`Failed to load logs: ${error.message}`)
      toast.error('Failed to load task logs')
    } finally {
      setLoadingLogs(false)
    }
  }

  const handleCopyLogs = () => {
    navigator.clipboard.writeText(taskLogs)
    toast.success('Logs copied to clipboard')
  }

  const columns: ColumnsType<Task> = [
    {
      title: 'Task Name',
      dataIndex: 'task_name',
      key: 'task_name',
      render: (name) => <span style={{ fontWeight: 500 }}>{name}</span>,
    },
    {
      title: 'Type',
      dataIndex: 'task_type',
      key: 'task_type',
      width: 120,
      render: (type) => type && <Tag color="blue">{type}</Tag>,
    },
    {
      title: 'Status',
      dataIndex: 'status',
      key: 'status',
      width: 120,
      render: (status: TaskStatus) => <StatusTag status={status} />,
    },
    {
      title: 'Progress',
      dataIndex: 'progress',
      key: 'progress',
      width: 150,
      render: (progress) => (
        <Progress percent={Math.round(progress)} size="small" status={progress === 100 ? 'success' : 'active'} />
      ),
    },
    {
      title: 'Created',
      dataIndex: 'created_at',
      key: 'created_at',
      width: 150,
      render: (date) => dayjs(date).format('YYYY-MM-DD HH:mm'),
    },
    {
      title: 'Actions',
      key: 'actions',
      width: 200,
      render: (_, record) => (
        <Space>
          <Tooltip title="View Logs">
            <Button
              type={record.status === 'failed' ? 'default' : 'text'}
              size="small"
              icon={<FileTextOutlined />}
              onClick={() => handleViewLogs(record)}
              danger={record.status === 'failed'}
            >
              Logs
            </Button>
          </Tooltip>
          {record.status === 'completed' && (
            <Tooltip title="View Results">
              <Button
                type="primary"
                size="small"
                icon={<EyeOutlined />}
                onClick={() => navigate(`/results/${record.id}`)}
              >
                Results
              </Button>
            </Tooltip>
          )}
          {record.status === 'running' && (
            <Tooltip title="Cancel Task">
              <Button size="small" danger icon={<StopOutlined />} onClick={() => handleCancelTask(record)}>
                Cancel
              </Button>
            </Tooltip>
          )}
        </Space>
      ),
    },
  ]

  const statisticItems: StatisticItem[] = [
    {
      key: 'total',
      title: 'Total',
      value: stats?.total_tasks || 0,
      valueStyle: { color: 'var(--color-primary)' },
      colSpan: { xs: 24, sm: 12, lg: 6 },
    },
    {
      key: 'running',
      title: 'Running',
      value: stats?.running_tasks || 0,
      valueStyle: { color: 'var(--color-warning)' },
      prefix: <SyncOutlined spin />,
      colSpan: { xs: 24, sm: 12, lg: 6 },
    },
    {
      key: 'completed',
      title: 'Completed',
      value: stats?.completed_tasks || 0,
      valueStyle: { color: 'var(--color-success)' },
      prefix: <CheckCircleOutlined />,
      colSpan: { xs: 24, sm: 12, lg: 6 },
    },
    {
      key: 'failed',
      title: 'Failed',
      value: stats?.failed_tasks || 0,
      valueStyle: { color: 'var(--color-error)' },
      prefix: <CloseCircleOutlined />,
      colSpan: { xs: 24, sm: 12, lg: 6 },
    },
  ]

  // Show skeleton on initial load
  if (initialLoad && loading) {
    return <PageSkeleton hasHeader rows={8} />
  }

  return (
    <div>
      {/* Statistics with staggered animation */}
      <StaggeredList staggerDelay={80} baseDelay={0} direction="up">
        <StatisticCard items={statisticItems} />
      </StaggeredList>

      {/* Header with filters */}
      <FadeIn direction="up" delay={100} duration={300}>
        <Space direction="vertical" size="middle" style={{ width: '100%', marginTop: 24, marginBottom: 16 }}>
          <Space align="center">
            <ThunderboltOutlined style={{ fontSize: 28, color: 'var(--color-primary)' }} />
            <div>
              <Title level={3} style={{ margin: 0 }}>
                Task Monitoring
              </Title>
              <Text type="secondary">Real-time pipeline execution tracking</Text>
            </div>
          </Space>

          <PageHeader
            left={
              <Select
                placeholder="All Projects"
                style={{ width: 300 }}
                value={selectedProject || undefined}
                onChange={setSelectedProject}
                allowClear
              >
                {items.map((p) => (
                  <Option key={p.id} value={p.id}>
                    {p.name}
                  </Option>
                ))}
              </Select>
            }
            right={
              <Button type="primary" icon={<PlusOutlined />}>
                New Task
              </Button>
            }
          />
        </Space>
      </FadeIn>

      {/* Task table with enhanced empty state */}
      <FadeIn direction="up" delay={200} duration={300}>
        {tasks.length === 0 && !loading ? (
          <EnhancedEmptyState
            type="noData"
            title="No tasks yet"
            description={
              selectedProject
                ? 'No tasks have been created for this project. Start a pipeline execution to create tasks.'
                : 'No tasks have been created yet. Select a project and execute a pipeline to get started.'
            }
            action={
              selectedProject
                ? {
                    text: 'Execute Pipeline',
                    onClick: () => navigate('/pipelines'),
                    icon: <ThunderboltOutlined />,
                  }
                : {
                    text: 'View Pipelines',
                    onClick: () => navigate('/pipelines'),
                    icon: <ThunderboltOutlined />,
                  }
            }
            size="default"
          />
        ) : (
          <DataTable
            columns={columns}
            dataSource={tasks}
            rowKey="id"
            loading={loading}
            pagination={{ pageSize: 20, showSizeChanger: true, showTotal: (total) => `Total ${total} tasks` }}
            emptyText="No Tasks"
            emptyDescription="No tasks have been created yet"
          />
        )}
      </FadeIn>

      {/* Task Logs Drawer */}
      <Drawer
        title={
          <Space>
            <FileTextOutlined />
            <span>Task Logs: {currentTask?.task_name}</span>
          </Space>
        }
        placement="right"
        width={720}
        open={logDrawerVisible}
        onClose={() => setLogDrawerVisible(false)}
        extra={
          <Button icon={<CopyOutlined />} onClick={handleCopyLogs} disabled={!taskLogs || loadingLogs}>
            Copy
          </Button>
        }
      >
        {/* Error message for failed tasks */}
        {currentTask?.status === 'failed' && currentTask?.error_message && (
          <Alert
            message="Task Failed"
            description={currentTask.error_message}
            type="error"
            showIcon
            icon={<CloseCircleOutlined />}
            style={{ marginBottom: 16 }}
          />
        )}

        {/* Task metadata */}
        <div style={{ marginBottom: 16, padding: 12, background: 'var(--color-gray-50)', borderRadius: 8 }}>
          <Space direction="vertical" size="small" style={{ width: '100%' }}>
            <div>
              <Text type="secondary">Task ID: </Text>
              <Text copyable={{ text: currentTask?.id || '' }} style={{ fontFamily: 'monospace', fontSize: 12 }}>
                {currentTask?.id}
              </Text>
            </div>
            <div>
              <Text type="secondary">Status: </Text>
              {currentTask && <StatusTag status={currentTask.status} />}
            </div>
            <div>
              <Text type="secondary">Progress: </Text>
              <Text strong>{Math.round(currentTask?.progress || 0)}%</Text>
            </div>
            {currentTask?.started_at && (
              <div>
                <Text type="secondary">Started: </Text>
                <Text>{dayjs(currentTask.started_at).format('YYYY-MM-DD HH:mm:ss')}</Text>
              </div>
            )}
            {currentTask?.completed_at && (
              <div>
                <Text type="secondary">Completed: </Text>
                <Text>{dayjs(currentTask.completed_at).format('YYYY-MM-DD HH:mm:ss')}</Text>
              </div>
            )}
          </Space>
        </div>

        {/* Log content */}
        <div
          style={{
            background: '#1e1e1e',
            color: '#d4d4d4',
            padding: 16,
            borderRadius: 8,
            fontFamily: 'monospace',
            fontSize: 13,
            lineHeight: 1.6,
            maxHeight: 'calc(100vh - 400px)',
            overflow: 'auto',
            whiteSpace: 'pre-wrap',
            wordBreak: 'break-word',
          }}
        >
          {loadingLogs ? (
            <div style={{ textAlign: 'center', padding: 40 }}>
              <Spin tip="Loading logs..." />
            </div>
          ) : (
            taskLogs || 'No logs available'
          )}
        </div>
      </Drawer>
    </div>
  )
}

export default TaskList
