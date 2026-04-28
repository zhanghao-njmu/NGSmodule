/**
 * Task List Page — TanStack Query + realtime migration.
 *
 * Replaces the previous Zustand store + manual WebSocket setup with:
 *   - `useTaskList` / `useTaskStatsSummary` (TanStack Query)
 *   - `useCancelTask` mutation with automatic cache invalidation
 *   - `useTaskLogs` (lazy, only fires when the drawer opens)
 *   - Realtime layer is already mounted globally in MainLayout, so list
 *     and detail caches refresh automatically when tasks progress.
 */
import type React from 'react'
import { useMemo, useState } from 'react'
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
import dayjs from 'dayjs'

import { useTaskList, useTaskStatsSummary, useCancelTask, useTaskLogs, useProjectList } from '@/hooks/queries'
import {
  PageHeader,
  DataTable,
  StatisticCard,
  StatusTag,
  PageSkeleton,
  FadeIn,
  EnhancedEmptyState,
  StaggeredList,
} from '@/components/common'
import { confirm } from '@/components/common/ConfirmDialog'
import { toast } from '@/utils/notification'
import type { StatisticItem } from '@/components/common'
import type { Task, TaskStatus } from '@/types/task'

const { Option } = Select
const { Title, Text } = Typography

export const TaskList: React.FC = () => {
  const navigate = useNavigate()

  const [selectedProject, setSelectedProject] = useState<string>('')
  const [logDrawerVisible, setLogDrawerVisible] = useState(false)
  const [currentTask, setCurrentTask] = useState<Task | null>(null)

  // ---- queries -------------------------------------------------------------

  // Project dropdown options
  const { data: projectList } = useProjectList({ limit: 200 })
  const projects = (projectList as any)?.items ?? (projectList as any)?.data ?? []

  // Task list — query key changes when filter changes, so refetch is automatic
  const taskListParams = selectedProject ? { project_id: selectedProject } : {}
  const { data: taskListData, isLoading: tasksLoading, isFetching: tasksFetching } = useTaskList(taskListParams)
  const tasks: Task[] = useMemo(() => (taskListData as any)?.items ?? (taskListData as any) ?? [], [taskListData])

  // Stats — kept fresh by realtime invalidation when tasks finish
  const { data: stats } = useTaskStatsSummary(selectedProject ? { project_id: selectedProject } : undefined)

  // Logs — only fetched while the drawer is open (enabled via taskId)
  const { data: logsResponse, isLoading: loadingLogs } = useTaskLogs(logDrawerVisible ? currentTask?.id : undefined)
  const taskLogs = (logsResponse as any)?.log_content || ''

  // ---- mutations -----------------------------------------------------------

  const cancelMutation = useCancelTask()

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
          await cancelMutation.mutateAsync(task.id)
          loadingToast()
          toast.success('任务已取消')
        } catch {
          loadingToast()
          toast.error('取消失败，请重试')
        }
      },
    })
  }

  const handleViewLogs = (task: Task) => {
    setCurrentTask(task)
    setLogDrawerVisible(true)
    // No manual fetch needed — useTaskLogs above is enabled by the drawer state
  }

  const handleCopyLogs = () => {
    navigator.clipboard.writeText(taskLogs)
    toast.success('Logs copied to clipboard')
  }

  // ---- table config --------------------------------------------------------

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
              <Button
                size="small"
                danger
                icon={<StopOutlined />}
                onClick={() => handleCancelTask(record)}
                loading={cancelMutation.isPending && cancelMutation.variables === record.id}
              >
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
      value: (stats as any)?.total_tasks || 0,
      valueStyle: { color: 'var(--color-primary)' },
      colSpan: { xs: 24, sm: 12, lg: 6 },
    },
    {
      key: 'running',
      title: 'Running',
      value: (stats as any)?.running_tasks || 0,
      valueStyle: { color: 'var(--color-warning)' },
      prefix: <SyncOutlined spin />,
      colSpan: { xs: 24, sm: 12, lg: 6 },
    },
    {
      key: 'completed',
      title: 'Completed',
      value: (stats as any)?.completed_tasks || 0,
      valueStyle: { color: 'var(--color-success)' },
      prefix: <CheckCircleOutlined />,
      colSpan: { xs: 24, sm: 12, lg: 6 },
    },
    {
      key: 'failed',
      title: 'Failed',
      value: (stats as any)?.failed_tasks || 0,
      valueStyle: { color: 'var(--color-error)' },
      prefix: <CloseCircleOutlined />,
      colSpan: { xs: 24, sm: 12, lg: 6 },
    },
  ]

  // Show skeleton on first load, but allow background refetches without skeleton
  if (tasksLoading && tasks.length === 0) {
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
                onChange={(v) => setSelectedProject(v ?? '')}
                allowClear
              >
                {projects.map((p: any) => (
                  <Option key={p.id} value={p.id}>
                    {p.name}
                  </Option>
                ))}
              </Select>
            }
            right={
              <Button type="primary" icon={<PlusOutlined />} onClick={() => navigate('/pipelines')}>
                New Task
              </Button>
            }
          />
        </Space>
      </FadeIn>

      {/* Task table */}
      <FadeIn direction="up" delay={200} duration={300}>
        {tasks.length === 0 && !tasksLoading ? (
          <EnhancedEmptyState
            type="noData"
            title="No tasks yet"
            description={
              selectedProject
                ? 'No tasks have been created for this project. Start a pipeline execution to create tasks.'
                : 'No tasks have been created yet. Select a project and execute a pipeline to get started.'
            }
            action={{
              text: selectedProject ? 'Execute Pipeline' : 'View Pipelines',
              onClick: () => navigate('/pipelines'),
              icon: <ThunderboltOutlined />,
            }}
            size="default"
          />
        ) : (
          <DataTable
            columns={columns}
            dataSource={tasks}
            rowKey="id"
            loading={tasksFetching && tasks.length === 0}
            pagination={{
              pageSize: 20,
              showSizeChanger: true,
              showTotal: (total) => `Total ${total} tasks`,
            }}
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
