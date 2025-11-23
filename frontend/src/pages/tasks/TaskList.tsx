/**
 * Task List Page - Task monitoring with real-time updates
 * Modernized with animations and enhanced UI components
 */
import React, { useEffect, useState } from 'react'
import { useNavigate } from 'react-router-dom'
import {
  Button,
  Tag,
  Progress,
  Space,
  Select,
  Tooltip,
  Typography,
} from 'antd'
import {
  PlusOutlined,
  SyncOutlined,
  CheckCircleOutlined,
  CloseCircleOutlined,
  StopOutlined,
  EyeOutlined,
  ThunderboltOutlined,
} from '@ant-design/icons'
import type { ColumnsType } from 'antd/es/table'
import { useTaskStore } from '../../store/taskStore'
import { useProjectStore } from '../../store/projectStore'
import { websocketService } from '../../services/websocket.service'
import {
  PageHeader,
  DataTable,
  StatisticCard,
  StatusTag,
  PageSkeleton,
  FadeIn,
  StaggeredList,
  EnhancedEmptyState,
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
  const {
    tasks,
    stats,
    loading,
    fetchTasks,
    fetchStats,
    cancelTask,
    handleWebSocketMessage,
  } = useTaskStore()
  const { items, fetchItems } = useProjectStore()
  const [selectedProject, setSelectedProject] = useState<string>('')
  const [initialLoad, setInitialLoad] = useState(true)

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
      }
    })
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
        <Progress
          percent={Math.round(progress)}
          size="small"
          status={progress === 100 ? 'success' : 'active'}
        />
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
      width: 150,
      render: (_, record) => (
        <Space>
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
    return <PageSkeleton hasHeader hasStats rows={8} />
  }

  return (
    <div>
      {/* Statistics with animation */}
      <FadeIn direction="up" delay={0} duration={300}>
        <StatisticCard items={statisticItems} />
      </FadeIn>

      {/* Header with filters */}
      <FadeIn direction="up" delay={100} duration={300}>
        <Space direction="vertical" size="middle" style={{ width: '100%', marginTop: 24, marginBottom: 16 }}>
          <Space align="center">
            <ThunderboltOutlined style={{ fontSize: 28, color: 'var(--color-primary)' }} />
            <div>
              <Title level={3} style={{ margin: 0 }}>
                Task Monitoring
              </Title>
              <Text type="secondary">
                Real-time pipeline execution tracking
              </Text>
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
    </div>
  )
}

export default TaskList
